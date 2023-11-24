
FFTW.set_num_threads(32)
theme(:dao)
ENV["GKSwstype"]="nul"
#ENV["LIBPOWSPEC_PATH"]="/home/users/d/dforeros/codes/libpowspec"
#using Powspec
ENV["OMP_NUM_THREADS"] = "32"
gr()
include("cosmo.jl")
include("bispec_wrap.jl")


function check_in_cell(x, y, z, ii, jj, kk, box_size, n_grid)
    bin_size = box_size / n_grid
    xmin = (ii - 1) * bin_size
    xmax = (ii ) * bin_size
    ymin = (jj - 1) * bin_size
    ymax = (jj ) * bin_size
    zmin = (kk - 1) * bin_size
    zmax = (kk ) * bin_size

    (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & (z >= zmin) & (z <= zmax)

end #func
@inline function check_in_range(x, ii, box_size, n_grid)
    bin_size = box_size / n_grid
    xmin = (ii - 1) * bin_size
    xmax = (ii ) * bin_size

    (x >= xmin) & (x <= xmax)

end #func
function check_in_cell_verbose(x, y, z, ii, jj, kk, box_size, n_grid)
    bin_size = typeof(x)(box_size / n_grid)
    xmin = (ii - 1) * bin_size
    xmax = (ii ) * bin_size
    ymin = (jj - 1) * bin_size
    ymax = (jj ) * bin_size
    zmin = (kk - 1) * bin_size
    zmax = (kk ) * bin_size
    @show xmin, x, xmax
    @show ymin, y, ymax
    @show zmin, z, zmax
    (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax) & (z >= zmin) & (z <= zmax)

end #func
function pos_to_grid(position, grid_size, box_size)
    1 + (Int(floor(grid_size * position / box_size) + grid_size) % grid_size )
end #func
function coordinate_separation(a, b, box_size)
    delta = abs(a - b)
    return (delta > 0.5*box_size ? delta - box_size : delta)*sign(a-b)
end


# Function that accumulates the histogram of distances
function build_histogram!(d2,hist, bin_edges)
    ibin = searchsortedlast(bin_edges, sqrt(d2))
    #ibin = floor(Int,d) + 1
    if (ibin > 0) & (ibin < length(bin_edges))
        hist[ibin] += 1
    end #if
    return hist
end;

function compute_2pcf(catalog, box_size, bin_edges)
    sides = [eltype(catalog)(b) for b in box_size]
    #cutoff = box_size[1] / size(counts_ref, 1)
    cutoff = eltype(catalog)(bin_edges[end])
    box = Box(sides, cutoff)
    
    cl = CellList(catalog, box)
    
    # Initialize (and preallocate) the histogram
    hist = zeros(Int,size(bin_edges,1)-1);
    println("Counting pairs...")
    # Run calculation
    CellListMap.map_pairwise!(
        (_,_,_,_,d2,hist) -> build_histogram!(d2,hist, bin_edges),
        hist,box,cl; show_progress = true, parallel = true
    )
    println("Done")
    #@show hist
    #hist = hist  / (length(satellites[1]) * (length(centrals[1])))
    N = size(catalog,2)
    hist = hist  / (N * (N - 1) / 2)
    #@show hist
    norm = 4 * pi .*  (bin_edges[2:end].^3 .- bin_edges[1:end-1].^3) ./ (box_size[1] * box_size[2] * box_size[3] * 3)
    #@show norm
    hist = hist ./ norm .- 1
    #@show hist
    centers = bin_edges[1:end-1] .+ 0.5 .* diff(bin_edges)
    centers, hist
end #

compute_2pcf(x, y, z, box_size, bin_edges) = compute_2pcf(Matrix(hcat(x,y,z)'), box_size, bin_edges)
    

function compute_2pcf(x1, y1, z1, x2, y2, z2, box_size, bin_edges)
    sides = [eltype(x1)(b) for b in box_size]
    #cutoff = box_size[1] / size(counts_ref, 1)
    cutoff = eltype(x1)(bin_edges[end])
    box = Box(sides, cutoff)
    
    cl = CellList(Matrix(hcat(x1,y1,z1)'), Matrix(hcat(x2,y2,z2)'), box)
        
    # Initialize (and preallocate) the histogram
    hist = zeros(Int,size(bin_edges,1)-1);    
    # Run calculation
    println("Counting pairs...")
    map_pairwise!(
        (x,y,i,j,d2,hist) -> build_histogram!(d2,hist, bin_edges),
        hist,box,cl
    )
    println("Done")
    hist = hist  / (length(x1[1]) * (length(x2[1])))
    norm = @. (4/3) * π * (bin_edges[2:end]^3 -bin_edges[1:end-1]^3) / (box_size[1] * box_size[2] * box_size[3])
    hist ./= norm
    centers = bin_edges[1:end-1] .+ 0.5 .* diff(bin_edges)
    centers, hist
end #func
@inline function wrap_indices(index, ngrid)
    if index > ngrid
        return index - ngrid
    elseif index < 1
        return index + ngrid
    else
        return index
    end #if
end #func
function ngp_assign(px::AbstractVector{T}, py::AbstractVector{T}, pz::AbstractVector{T}, boxsize, ngrid) where T<:Real

    counts = zeros(UInt32, (ngrid,ngrid,ngrid))
    ind3d = zeros(UInt32, size(px))
    aa = zero(ind3d)
    bb = zero(ind3d)
    cc = zero(ind3d)
    
    for ii in eachindex(px)
        i = wrap_indices(1 + Int(floor(ngrid * px[ii] / boxsize)), ngrid)
        j = wrap_indices(1 + Int(floor(ngrid * py[ii] / boxsize)), ngrid)
        k = wrap_indices(1 + Int(floor(ngrid * pz[ii] / boxsize)), ngrid)
        counts[i,j,k] += 1
        ind3d[ii] = LinearIndices(counts)[i,j,k]
        aa[ii] = i
        bb[ii] = j
        cc[ii] = k
    end #for
    return ind3d, counts, aa, bb, cc
end #func

function read_dm_particles(box_size)
    dm_particles = [Array{Float32}(undef, 360^3) for _ in 1:3]
    grid_particles = [Array{Float32}(undef, 360^3) for _ in 1:3]
    buffer = zero(dm_particles[1])
    axes = ['x', 'y', 'z']
    for i = 1:3
        read!("/home/users/d/dforeros/projects/desi-patchy/run/pos$(axes[i]).dat", dm_particles[i])
        grid_particles[i] .= dm_particles[i]
        read!("/srv/beegfs/scratch/users/d/dforeros/projects/desi-patchy/run/BOXpos$(axes[i])OM0.315OL0.685G360V2000.0_ALPTrs6.000z0.500.dat", buffer)
        #read!("/home/users/d/dforeros/projects/desi-patchy/run/PSI$(axes[i])z0.500.dat", buffer)
        @. dm_particles[i] = (buffer + box_size) % box_size
    end #for
    dm_particles, grid_particles
end #func

function read_displacement()
    displacements = [Array{Float32}(undef, 360^3) for _ in 1:3]
    axes = ['x', 'y', 'z']
    for i = 1:3
        read!("/home/users/d/dforeros/projects/desi-patchy/run/PSI$(axes[i])z0.500.dat", displacements[i])
    end #for
    displacements
end #func

function read_velocity()
    displacements = [Array{Float32}(undef, 360^3) for _ in 1:3]
    axes = ['x', 'y', 'z']
    for i = 1:3
        read!("/home/users/d/dforeros/projects/desi-patchy/run/VE$(axes[i])z0.500.dat", displacements[i])
    end #for
    displacements
end #func

function read_eulerian_velocity()
    displacements = [Array{Float32}(undef, 360^3) for _ in 1:3]
    axes = ['x', 'y', 'z']
    for i = 1:3
        read!("/home/users/d/dforeros/projects/desi-patchy/run/VE$(axes[i])EULz0.500.dat", displacements[i])
    end #for
    displacements
end #func

function read_disp_zeld()
    displacements = [Array{Float32}(undef, 360^3) for _ in 1:3]
    axes = ['x', 'y', 'z']
    for i = 1:3
        read!("/home/users/d/dforeros/projects/desi-patchy/run/PSIZELD$(axes[i]).dat", displacements[i])
    end #for
    displacements
end #func

function read_cosmic_web()
    buffer = Array{Float32}(undef, 360^3)
    #read!("/home/users/d/dforeros/projects/desi-patchy/run/TwebDelta_OM0.315OL0.685G360V2000.0z0.500.dat", buffer)
    read!("/srv/beegfs/scratch/users/d/dforeros/projects/desi-patchy/run/TwebDelta_OM0.315OL0.685G360V2000.0z0.500.dat", buffer)
    buffer
end #func

function read_dm_dens()
    buffer = Array{Float32}(undef, 360^3)
    #read!("/home/users/d/dforeros/projects/desi-patchy/run/TwebDelta_OM0.315OL0.685G360V2000.0z0.500.dat", buffer)
    read!("/srv/beegfs/scratch/users/d/dforeros/projects/desi-patchy/run/deltaBOXOM0.315OL0.685G360V2000.0_ALPTrs6.000z0.500.dat", buffer)
    buffer
end #func




function read_h5(filename, columns)
    store = h5open(filename, "r")
    @show keys(store)
    catalog = DataFrame(reduce(hcat, [map(Float32, store[key][:]) for key in columns]), map(lowercase, columns))
    close(store)
    println("Done reading H5")
    catalog
end #fun


function apply_rsd!(new_pos, pos, vel, redshift, cosmo::Cosmology)
    efunc = E(cosmo, redshift)
    Threads.@threads for i = eachindex(vel)
        disp = vel[i] * eltype(pos)((1 + redshift) / (100. * efunc))
        #println(vel[i])
        new_pos[i] = pos[i] + disp
    end #func
    new_pos
end #func


function read_cic!(output::AbstractVector{T}, field::AbstractArray{T, 3}, pos::AbstractVector{<:AbstractVector{T}}, box_size::AbstractVector, box_min::AbstractVector; wrap = true )  where T <: Real

    dims = size(field)
    cell_size = map(T, box_size ./ dims)
    
    @inbounds Threads.@threads for i in eachindex(pos[1])
        dist_x = (pos[1][i] - box_min[1]) / cell_size[1]
        dist_y = (pos[2][i] - box_min[2]) / cell_size[2]
        dist_z = (pos[3][i] - box_min[3]) / cell_size[3]

        dist_i_x = Int(floor(dist_x)) 
        dist_i_y = Int(floor(dist_y)) 
        dist_i_z = Int(floor(dist_z)) 

        ux = dist_x - dist_i_x
        uy = dist_y - dist_i_y
        uz = dist_z - dist_i_z

        dx = 1 - ux
        dy = 1 - uy
        dz = 1 - uz

        dist_i_x += 1
        dist_i_y += 1
        dist_i_z += 1


        index_d_x = (dist_i_x > dims[1]) & wrap ? dist_i_x - dims[1] : dist_i_x
        index_d_y = (dist_i_y > dims[2]) & wrap ? dist_i_y - dims[2] : dist_i_y
        index_d_z = (dist_i_z > dims[3]) & wrap ? dist_i_z - dims[3] : dist_i_z

        index_u_x = index_d_x + 1
        index_u_y = index_d_y + 1
        index_u_z = index_d_z + 1


        index_u_x = (index_u_x > dims[1]) & wrap ? index_u_x - dims[1] : index_u_x
        index_u_y = (index_u_y > dims[2]) & wrap ? index_u_y - dims[2] : index_u_y
        index_u_z = (index_u_z > dims[3]) & wrap ? index_u_z - dims[3] : index_u_z


        output[i] = field[index_d_x, index_d_y, index_d_z] * dx * dy * dz + 
                    field[index_d_x, index_d_y, index_u_z] * dx * dy * uz +
                    field[index_d_x, index_u_y, index_d_z] * dx * uy * dz +
                    field[index_d_x, index_u_y, index_u_z] * dx * uy * uz +
                    field[index_u_x, index_d_y, index_d_z] * ux * dy * dz +
                    field[index_u_x, index_d_y, index_u_z] * ux * dy * uz +
                    field[index_u_x, index_u_y, index_d_z] * ux * uy * dz +
                    field[index_u_x, index_u_y, index_u_z] * ux * uy * uz

    end #for
    output
end #func
read_cic!(output::AbstractVector{T}, field::AbstractArray{T, 3}, catalog::AbstractArray{T,2}, box_size::AbstractVector, box_min::AbstractVector; wrap = true )  where T <: Real = read_cic!(output, field, [@view(catalog[i,:]) for i=1:3], box_size, box_min; wrap = wrap )
function read_cic(field::AbstractArray{T,3}, position::AbstractVector{T}, box_size::AbstractVector, box_min::AbstractVector; wrap = true) where T<: Real
    dims = size(field)
    cell_size = map(T, box_size ./ dims)
    dist_x = 1 + ((position[1] - box_min[1]) / cell_size[1])
    dist_y = 1 + ((position[2] - box_min[2]) / cell_size[2])
    dist_z = 1 + ((position[3] - box_min[3]) / cell_size[3])

    dist_i_x = Int(floor(dist_x)) 
    dist_i_y = Int(floor(dist_y)) 
    dist_i_z = Int(floor(dist_z)) 

    ux = dist_x - dist_i_x
    uy = dist_y - dist_i_y
    uz = dist_z - dist_i_z

    dx = 1 - ux
    dy = 1 - uy
    dz = 1 - uz

    index_d_x = wrap ? wrap_indices(dist_i_x, dims[1]) : dist_i_x
    index_d_y = wrap ? wrap_indices(dist_i_y, dims[2]) : dist_i_y
    index_d_z = wrap ? wrap_indices(dist_i_z, dims[3]) : dist_i_z
    
    index_u_x = index_d_x + 1
    index_u_y = index_d_y + 1
    index_u_z = index_d_z + 1

    index_u_x = wrap ? wrap_indices(index_u_x, dims[1]) : index_u_x
    index_u_y = wrap ? wrap_indices(index_u_y, dims[2]) : index_u_y
    index_u_z = wrap ? wrap_indices(index_u_z, dims[3]) : index_u_z

    
    

    output = field[index_d_x, index_d_y, index_d_z] * dx * dy * dz + 
             field[index_d_x, index_d_y, index_u_z] * dx * dy * uz +
             field[index_d_x, index_u_y, index_d_z] * dx * uy * dz +
             field[index_d_x, index_u_y, index_u_z] * dx * uy * uz +
             field[index_u_x, index_d_y, index_d_z] * ux * dy * dz +
             field[index_u_x, index_d_y, index_u_z] * ux * dy * uz +
             field[index_u_x, index_u_y, index_d_z] * ux * uy * dz +
             field[index_u_x, index_u_y, index_u_z] * ux * uy * uz
    output
end #func
function renormalize_alpt_vel!(vels, cosmo, z)
    println("Renormalizing velocities")
    hconst  = 1  #Will be cancelled
    cvel = 1. / (cgs_km / cgs_Mpc)
    #H_eval = H(cosmo, z)# * cgs_km / cgs_Mpc / cgs_sec
    H_eval = 100. * E(cosmo, z)
    #@show H_eval
    cpvel = growth_rate_approx(cosmo, z) * (H_eval ) / (1 + z)
    #D1 = growth_factor_approx(cosmo, z) / growth_factor_approx(cosmo, 0)
    D1 = 1
    Threads.@threads for i = eachindex(vels[1])
        for j =1:3
            vels[j][i] *= D1 * cvel / hconst * cpvel
            vels[j][i] *= cgs_km / cgs_sec * hconst
        end #for
    end  #for
    vels
end #func



function k_vec(field::AbstractArray{T}, box_size::AbstractVector) where T<:Real
    dims = [size(field)...]
    sample_rate = map(T, 2π .* dims ./ box_size)
    kx = rfftfreq(dims[1], sample_rate[1])
    ky = fftfreq(dims[2], sample_rate[2])
    kz = fftfreq(dims[3], sample_rate[3])
    (kx, ky, kz)
end #func

function σ(x::T) where T
    one(T) / (one(T) + exp(-x))
end #func

function smooth!(field::AbstractArray{T, 3}, smoothing_radius::T, box_size::AbstractVector; fft_plan = nothing) where T <: Real

    if fft_plan == nothing
        fft_plan = plan_rfft(field)
    end #if
    field_k = fft_plan * field
    k⃗ = k_vec(field, box_size)
    @inbounds Threads.@threads for I in CartesianIndices(field_k)
        k² = k⃗[1][I[1]]^2 + k⃗[2][I[2]]^2 + k⃗[3][I[3]]^2
        #field_k[I] *= exp(-0.5 * smoothing_radius^2 * k²)
        #field_k[I] *= exp(0.5 * smoothing_radius^2 * k²)
        field_k[I] *= 1 + σ(0.6) * smoothing_radius * k²
    end #for
    ldiv!(field, fft_plan, field_k)
    field
end #func

function new_smooth!(field::AbstractArray{T, 3}, smoothing_radius::T, box_size::AbstractVector; fft_plan = nothing) where T <: Real

    if fft_plan == nothing
        fft_plan = plan_rfft(field)
    end #if
    field_k = fft_plan * field
    k⃗ = k_vec(field, box_size)
    @inbounds Threads.@threads for I in CartesianIndices(field_k)
        k² = k⃗[1][I[1]]^2 + k⃗[2][I[2]]^2 + k⃗[3][I[3]]^2
        field_k[I] *= (1 + smoothing_radius * k²)
        #field_k[I] *= exp(0.5 * smoothing_radius^2 * k²)
        
    end #for
    ldiv!(field, fft_plan, field_k)
    field
end #func