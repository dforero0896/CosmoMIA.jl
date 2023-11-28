

function c_to_julia_order(arr, dims)
    permutedims(reshape(arr, (dims, dims, dims)), (3,2,1))
end #func

mutable struct ALPTResult{T}
    box_size::AbstractVector
    pos::AbstractVector{<:AbstractVector{T}}
    vel_field::AbstractVector{<:AbstractArray{T,3}}
    disp_field::AbstractVector{<:AbstractArray{T,3}} 
    delta_field::AbstractArray{T,3} 
    dweb_field::AbstractArray{UInt16,3}
    function ALPTResult{T}( box_size, ngrid; 
                        pos_path = nothing, 
                        vel_path = nothing, 
                        disp_path = nothing, 
                        delta_path = nothing,
                        dweb_path = nothing) where T<:Real
        if pos_path != nothing
            pos = [Array{T}(undef, ngrid^3) for _=1:3]
        else 
            pos = [Array{T}(undef, 0)]
        end #if
        if vel_path != nothing
            vel_field = [Array{T}(undef, ngrid^3) for _=1:3]
        else
            vel_field = [Array{T,3}(undef, 0,0,0)]
        end #if
        if disp_path != nothing
            disp_field = [Array{T}(undef, ngrid^3) for _=1:3]
        else 
            disp_field = [Array{T,3}(undef, 0,0,0)]
        end
        if delta_path != nothing
            delta_field = Array{T}(undef, ngrid^3)
            read!(delta_path, delta_field)
            delta_field = c_to_julia_order(delta_field, ngrid)
        else
            dweb_field = Array{T,3}(undef, 0,0,0)
        end
        if dweb_path != nothing
            dweb_field = Array{T}(undef, ngrid^3)
            read!(dweb_path, dweb_field)
            dweb_field = UInt16.(c_to_julia_order(dweb_field, ngrid))
        else
            dweb_field = Array{UInt16,3}(undef, 0,0,0)
        end
        axes = ['x', 'y', 'z']
        for i = 1:3
            if pos_path != nothing
                read!(Printf.format(Printf.Format(pos_path), axes[i]), pos[i])
                pos[i] .= (pos[i] .+ box_size[i]) .% box_size[i]
            end 
            if vel_path != nothing
                read!(Printf.format(Printf.Format(vel_path), axes[i]), vel_field[i])
            end 
            if disp_path != nothing
                read!(Printf.format(Printf.Format(disp_path), axes[i]), disp_field[i])
            end 
        end #for
        if vel_path != nothing
            vel_field = c_to_julia_order.(vel_field, ngrid)
        end 
        if disp_path != nothing
            disp_field = c_to_julia_order.(disp_field, ngrid)
        end         
        new{T}(box_size, pos, vel_field, disp_field, delta_field, dweb_field)
    end #func
end #struct

mutable struct SubgridCatalog{T}
    pos::AbstractMatrix{T}
    vel::AbstractMatrix{T}
    is_dm::BitVector
    collapse_to_idx::AbstractVector{UInt32}
    dweb::AbstractVector{UInt16}
    δ_dm::AbstractVector{T}
    r_min::AbstractVector{T}
    δ_max::AbstractVector{T}
    is_attractor::BitVector
    γ_par
    collapse_frac
    collapse_radius
end #struct

function copy(cat::SubgridCatalog{T}) where T
    println("Copying SubgridCatalog")
    SubgridCatalog(map(deepcopy, (cat.pos, cat.vel, 
                                  cat.is_dm,
                                  cat.collapse_to_idx, 
                                  cat.dweb, 
                                  cat.δ_dm, cat.r_min, 
                                  cat.δ_max, cat.is_attractor,
                                  cat.γ_par,
                                  cat.collapse_frac,
                                  cat.collapse_radius))...)
end #func
function copy!(target_cat::CT, source_cat::CT) where {CT, T}
    println("Copying SubgridCatalog in place")
    target_cat.pos .= source_cat.pos
    target_cat.vel .= source_cat.vel
    target_cat.is_dm .= source_cat.is_dm
    target_cat.collapse_to_idx .= source_cat.collapse_to_idx
    target_cat.dweb .= source_cat.dweb
    target_cat.δ_dm .= source_cat.δ_dm
    target_cat.r_min = source_cat.r_min
    target_cat.δ_max = source_cat.δ_max
    target_cat.is_attractor .= source_cat.is_attractor
    target_cat.γ_par = source_cat.γ_par
    target_cat.collapse_frac = source_cat.collapse_frac
    target_cat.collapse_radius = source_cat.collapse_radius
    target_cat
end #func
function assign_particles_to_gals(dm_particles, target_ncount::AbstractArray{IT, 3}, box_size, box_min, dm_cw_type, dm_dens, displacement, velocities, dist; debug=false) where IT <: Unsigned
    println("Assigning DM to galaxies...")
    grid_size = size(target_ncount)
    bin_size = box_size ./ grid_size
    dm_per_cell = [Vector{UInt32}() for _ = 1:prod(grid_size)]
    LI = LinearIndices(target_ncount)
    if debug
        number_dm = zero(target_ncount)
    end #if
    @inbounds for ii = eachindex(dm_particles[1])
        i = wrap_indices(1 + Int(floor(grid_size[1] * dm_particles[1][ii] / box_size[1])), grid_size[1])
        j = wrap_indices(1 + Int(floor(grid_size[2] * dm_particles[2][ii] / box_size[2])), grid_size[2])
        k = wrap_indices(1 + Int(floor(grid_size[3] * dm_particles[3][ii] / box_size[3])), grid_size[3])
        index_3d = LI[i, j, k]
        push!(dm_per_cell[index_3d], ii)
        if debug
            number_dm[i, j, k] += 1
        end #if
    end #for
    pos = Vector{SVector{3,eltype(dm_particles[1])}}()
    vel = Vector{SVector{3,eltype(velocities[1])}}()
    is_dm = BitVector()
    dweb = Vector{eltype(dm_cw_type)}()
    δ_dm = Vector{eltype(dm_dens)}()
    used_dm = 0
    sampled_new = 0

    for I = CartesianIndices(target_ncount)
        number_cen_in_cell = target_ncount[I]
        if number_cen_in_cell === zero(typeof(number_cen_in_cell))
            continue
        end #if
        TI = Tuple(I)
        grid_center = [(i - 0.5) * bin_size[1] for i = TI]
        L = LI[I]
        L_C = ((TI[3] - 1) + grid_size[3] * ((TI[2] - 1) + grid_size[2] * (TI[1] - 1))) + 1
        dm_particles_in_cell = dm_per_cell[L]
        number_dm_in_cell = length(dm_particles_in_cell)
        if debug
            @assert number_dm_in_cell == number_dm[I]
        end #if

        missing_counter = 0
        for par = 1:number_cen_in_cell
            if par <= number_dm_in_cell
                push!(pos, @SVector [dm_particles[ax][dm_particles_in_cell[par]] for ax = 1:3])
                push!(is_dm, true)
                used_dm += 1
            else
                missing_counter = par - number_dm_in_cell
                new_pos = Vector{eltype(dm_particles[1])}(undef, 3)
                lag_pos = Vector{eltype(dm_particles[1])}(undef, 3)
                if number_dm_in_cell > 0
                    missing_counter = missing_counter > number_dm_in_cell ? rand(1:number_dm_in_cell) : missing_counter
                    for ax = 1:3
                        lag_pos[ax] = (missing_counter <= number_dm_in_cell) ? dm_particles[ax][dm_particles_in_cell[missing_counter]] : grid_center[ax]
                        lag_pos[ax] -= displacement[ax][dm_particles_in_cell[missing_counter]]
                        lag_pos[ax] += rand(dist)
                    end #for
                    for ax = 1:3
                        #new_disp = read_cic(displacement[ax], lag_pos, box_size, box_min; wrap = true)
                        new_disp = displacement[ax][dm_particles_in_cell[missing_counter]]
                        new_pos[ax] = (lag_pos[ax] + new_disp + box_size[ax]) % box_size[ax]
                    end #for
                    
                else
                    for ax = 1:3
                        draw = rand(Uniform(-1, 1))
                        new_pos[ax] = (grid_center[ax] + (0.5 * bin_size[ax] * sign(draw) * (1 - sqrt(abs(draw)))) + box_size[ax]) % box_size[ax]   
                    end #for
                    sampled_new += 1
                end #if
                push!(pos, SVector{3}(new_pos))
                sampled_new += 1
                push!(is_dm, false)
            end #if
            push!(dweb, dm_cw_type[I])
            push!(δ_dm, read_cic(dm_dens, pos[end], box_size, box_min; wrap = true))
            #push!(vel, @SVector [velocities[ax][I] for ax = 1:3])
            push!(vel, @SVector [read_cic(velocities[ax], pos[end], box_size, box_min; wrap = true) for ax = 1:3])
        end #for
    end #for
    println("Used ", used_dm, " DM particles")
    println("Sampled ", sampled_new, " particles")
    println("Sampled ", round(100 * sampled_new / (used_dm + sampled_new)), "% of target_ncount")
    r_min = fill!(similar(δ_dm), +Inf)
    δ_max = fill!(similar(δ_dm), -Inf)
    collapse_to_idx = fill!(Vector{UInt32}(undef, size(δ_max,1)), 0)
    SubgridCatalog(reduce(hcat, pos), reduce(hcat,vel), is_dm, collapse_to_idx, dweb, δ_dm, r_min, δ_max, is_attractor_fun.(is_dm, dweb), 0f0, 1f0, 15f0)
    
end #func


assign_particles_to_gals(alpt_storage::ALPTResult{T}, 
                         target_ncount::AbstractArray{IT, 3}, 
                         box_size, box_min, dist; debug=false
                         ) where {T<:Real, IT<:Unsigned} = assign_particles_to_gals(alpt_storage.pos, 
                                                              target_ncount, 
                                                              box_size, 
                                                              box_min, 
                                                              alpt_storage.dweb_field, 
                                                              alpt_storage.delta_field, 
                                                              alpt_storage.disp_field, 
                                                              alpt_storage.vel_field, 
                                                              dist; debug=debug)



sinx(cos2x) = cos2x ≈ 1 ? 0 : sqrt(1 - cos2x)
collapse_func(r, par) = r * par
function collapse!(sat_pos, cen_pos, sat_id, cen_id, d2, catalog::SubgridCatalog{T}) where T
    collapse_frac = catalog.collapse_frac
    γ = catalog.γ_par
    r = sqrt(d2)
    cosϕ = (sat_pos[3] - cen_pos[3]) / r
    sinϕ = sinx(cosϕ^2)
    cosθ = sinϕ == 0 ? √2 / 2 : (sat_pos[1] - cen_pos[1])/ (r * sinϕ)
    sinθ = sinx(cosθ^2)
    r = collapse_func(r, collapse_frac)
    
    catalog.pos[1, sat_id] = cen_pos[1] + r * cosθ * sinϕ
    catalog.pos[2, sat_id] = cen_pos[2] + r * sinθ * sinϕ
    catalog.pos[3, sat_id] = cen_pos[3] + r * cosϕ

    
    δ = catalog.δ_dm[cen_id]
    for ax = 1:3
        draw = rand(Normal(0f0,1f0))
        catalog.vel[ax, sat_id] = catalog.vel[ax, sat_id] + 
                                    draw * γ * 1f1 * 
                                    ((1 + ifelse(δ <= 0, 0, δ)))^0.5 
    end #for
    catalog
end #func

is_attractor_fun(is_dm, cw_type) = (cw_type < 4) && is_dm

function inner_collapse_dm_dm!(cen_pos, sat_pos, cen_id, sat_id, d2, catalog::SubgridCatalog{T}) where T<:Real
    if sqrt(d2) > catalog.collapse_radius 
        return catalog
    end #if
    if ((catalog.is_attractor[cen_id] &&
        (catalog.is_dm[sat_id])) && 
        (sat_id != cen_id) && 
        (d2 > 1e-4))
        if d2 < catalog.r_min[sat_id]
            catalog = collapse!(sat_pos, cen_pos, sat_id, cen_id, d2, catalog::SubgridCatalog{T})
            #catalog.collapse_to_idx[sat_id] = UInt32(cen_id)
            catalog.r_min[sat_id] = d2
        end #if
    elseif ((catalog.is_attractor[sat_id] &&
            (catalog.is_dm[cen_id])) && 
            (sat_id != cen_id) && 
            (d2 > 1e-4))
            if d2 < catalog.r_min[cen_id]
                catalog = collapse!(cen_pos, sat_pos, cen_id, sat_id, d2, catalog::SubgridCatalog{T})
                #catalog.collapse_to_idx[sat_id] = UInt32(cen_id)
                catalog.r_min[cen_id] = d2
            end #if
    end #if
    catalog
end #func


function inner_collapse_ran_dm!(cen_pos, sat_pos, cen_id, sat_id, d2, catalog::SubgridCatalog{T}) where T<:Real
    if sqrt(d2) > catalog.collapse_radius 
        return catalog
    end #if
    if ((catalog.is_attractor[cen_id] && 
        !catalog.is_attractor[sat_id]) && 
        (d2 > 1e-4))
        if d2 < catalog.r_min[sat_id] || catalog.δ_max[sat_id] < catalog.δ_dm[cen_id]
            catalog = collapse!(sat_pos, cen_pos, sat_id, cen_id, d2, catalog)
            #catalog.collapse_to_idx[sat_id] = UInt32(cen_id)
            catalog.r_min[sat_id] = d2
            catalog.δ_max[sat_id] = catalog.δ_dm[cen_id]
        end #if
    elseif ((catalog.is_attractor[sat_id] && 
            !catalog.is_attractor[cen_id]) && 
            (d2 > 1e-4))
            if d2 < catalog.r_min[cen_id] || catalog.δ_max[cen_id] < catalog.δ_dm[sat_id]
                catalog = collapse!(cen_pos, sat_pos, cen_id, sat_id, d2, catalog)
                #catalog.collapse_to_idx[sat_id] = UInt32(cen_id)
                catalog.r_min[cen_id] = d2
                catalog.δ_max[cen_id] = catalog.δ_dm[sat_id]
            end #if
            
    end #if
    catalog
end #func

function reduce_dm_dm(output,output_threaded)
    catalog = output_threaded[1]
    for i in 2:length(output_threaded)
        if output_threaded[i].r_min < mind[3]
            mind = output_threaded[i]
        end
    end
    return mind
end


function subgrid_collapse(catalog::SubgridCatalog{T},
        params::AbstractVector,
        box_size::AbstractVector;
        show_progress = true) where T
    println("Starting collapse...")
    #catalog = copy_output(catalog)
    box = Box(box_size, maximum([8., params[1]]))
    cl = CellList(catalog.pos, box)
    catalog.collapse_radius = params[1]
    catalog.collapse_frac = params[3]
    catalog.γ_par = params[end]
    CellListMap.map_pairwise!((x, y, i, j, d2, catalog) -> inner_collapse_dm_dm!(x, y, i, j, d2, catalog),
                              catalog, box, cl, parallel=false, show_progress=show_progress)
    println("First collapse done")
    box = Box(box_size, maximum([8., params[2]]))
    cl = CellList(catalog.pos, box)
    catalog.collapse_radius = params[2]
    catalog.collapse_frac = params[4]
    CellListMap.map_pairwise!((x, y, i, j, d2, catalog) -> inner_collapse_ran_dm!(x, y, i, j, d2, catalog),
                               catalog, box, cl, parallel=false, show_progress=show_progress)

    catalog

end #func