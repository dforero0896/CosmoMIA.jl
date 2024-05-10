using Revise
using CosmoMIA
using CosmoCorr
using Plots
using Plots.PlotMeasures
using Printf
using FITSIO
using DataFrames
using StaticArrays
using Random
using Statistics
using Distributions
using StatsBase
using NPZ
using LaTeXStrings
using Optim
using FFTW
using CellListMap
FFTW.set_num_threads(32)
theme(:dao)
ENV["GKSwstype"]="nul"
ENV["OMP_NUM_THREADS"] = "32"
gr()
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--pos_path"
        default = "/home/astro/dforero/projects/desi-patchy/run-ridged/BOXpos%sOM0.315OL0.685G360V2000.0_ALPTrs6.000z1.100.dat"
        required = false
        "--vel_path"
        default = "/home/astro/dforero/projects/desi-patchy/run-ridged/VE%sEULz1.100.dat"
        required = false
        "--disp_path"
        default = "/home/astro/dforero/projects/desi-patchy/run-ridged/PSI%sz1.100.dat"
        required = false
        "--delta_path"
        default = "/home/astro/dforero/projects/desi-patchy/run-ridged/deltaBOXOM0.315OL0.685G360V2000.0_ALPTrs6.000z1.100.dat"
        required = false
        "--cweb_path"
        default = "/home/astro/dforero/projects/desi-patchy/run-ridged/TwebDelta_OM0.315OL0.685G360V2000.0lthD0.050z1.100.dat"
        required = false
        "--ref_catalog"
        help = "Reference catalog if fitting is to be made"
        required = false
        default = "/srv/astro/projects/cosmo3d/desi/SecondGenMocks/AbacusHOD/ELG/z1.100/AbacusSummit_base_c000_ph000/ELG_real_space.fits"
        "--do_fit"
        action = :store_true
        "--plot"
        action = :store_true
        "--skip_collapse"
        action = :store_true
        "--redshift"
        default = 1.1
        "--box_size"
        default = 2000.
        "--vel_boost"
        default = 8f0
        "--disp_smooth"
        default = 0.7f0
        "--step_size"
        default = -10f0
        "--number_counts"
        default = "/home/astro/dforero/cosmo3d/desi/patchy/ELG/z1.100/calibration/n_mock.dat"
        "--params"
        #default = [5.38945528422064, 5.343785994677935, 0.0003635007999382034, 0.3902633911503999, 0.3761037446458434]
        default = [6.7, 1., 0.4, 0.1, 0.2]
        arg_type = Float64
        nargs = 5
        "--output_path"
        required = false
        default = "./mock_out.npy"
    end #arg_table

    parsed_args = parse_args(s)
end #func


function fit_reference(args)
    cosmo = CosmoMIA.DESICosmology()
    box_size = @SVector [args["box_size"] for _ in 1:3]
    box_min = @SVector [0.0f0 for _ in 1:3]
    grid_size = Tuple(360 for _ in 1:3)
    k_edges = range(3f-3, 0.505, step = 5f-3)
    s_edges = Float32.(10 .^range(-1, log10(20), 51))
    los = [0,0,1]
    ref_size = 360
    bin_size = box_size ./ ref_size
    Random.seed!(42)

    println("Reading reference catalog")
    catalog = Float32.(DataFrame(FITS(args["ref_catalog"])[2]))
    println("Applying RSD to ref.")
    catalog[!,:zrsd] = CosmoMIA.apply_rsd!(zeros(eltype(catalog[!,:z]), size(catalog[!,:z],1)), catalog[!,:z], catalog[!,:vz], args["redshift"], cosmo)
    catalog[!, [:x, :y, :z, :zrsd]] .+= Float32(box_size[1] / 2)
    catalog[!, [:x, :y, :z, :zrsd]] .+= Float32(box_size[1])
    catalog[!, [:x, :y, :z, :zrsd]] .%= Float32(box_size[1])
    catalog[!,:w] = fill(1f0, size(catalog,1))

    println("Computing reference clustering")
    k, power_ref, _ = power_spectrum([catalog[!,:x], catalog[!,:y], catalog[!,:z]], catalog[!,:w], [0,2,4], k_edges, los, grid_size, box_size, box_min, 2)   
    k, power_ref_rsd, _ = power_spectrum([catalog[!,:x], catalog[!,:y], catalog[!,:zrsd]], catalog[!,:w], [0,2,4], k_edges, los, grid_size, box_size, box_min, 2)
    s, xi_ref = CosmoMIA.compute_2pcf(catalog[!,:x], catalog[!,:y], catalog[!,:zrsd], box_size, s_edges)
    if args["plot"]
        p = plot(k, k .^ 1.0 .* real.(@view(power_ref[1, :])), xscale=:identity, ylabel=L"P(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Ref.")
        p0 = plot(k, k .^ 1.0 .* real.(@view(power_ref_rsd[1, :])), xscale=:identity, ylabel=L"P_z(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Ref. l=0")
        p0 = plot!(p0, k, k .^ 1.0 .* real.(@view(power_ref_rsd[2, :])), xscale=:identity, ylabel=L"P_z(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Ref. l=2")
        p0 = plot!(p0, k, k .^ 1.0 .* real.(@view(power_ref_rsd[3, :])), xscale=:identity, ylabel=L"P_z(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Ref. l=4")
        px = plot(s, s .^ 0 .* xi_ref, xscale = :log10, yscale = :log10, label = "Ref.", xlabel = L"s~[{\rm  Mpc}/h]", ylabel = L"\xi(s)")
    end #if
    println("Creating reference field")
    ref_field = zeros(Float32, ref_size, ref_size, ref_size)
    CosmoCorr.ngp!(ref_field, catalog[!,:x], catalog[!,:y], catalog[!,:z], catalog[!,:w], box_size, box_min; wrap = true)
    ref_field = UInt32.(round.(ref_field))
    catalog = nothing

    println("Loading ALPT results")
    alpt_storage = CosmoMIA.ALPTResult{Float32}(box_size, ref_size; 
                                       pos_path = args["pos_path"], 
                                       vel_path = args["vel_path"], 
                                       disp_path = args["disp_path"], 
                                       delta_path = args["delta_path"],
                                       dweb_path = args["cweb_path"])
    println("Pre-boosting velocities.")
    CosmoMIA.renormalize_alpt_vel!(alpt_storage.vel_field, cosmo, args["redshift"])
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[1], args["vel_boost"], box_size)
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[2], args["vel_boost"], box_size)
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[3], args["vel_boost"], box_size)
    
    

    patchy_field = permutedims(reshape(read!(args["number_counts"], Array{Float32}(undef, ref_size^3)), (ref_size, ref_size, ref_size)), (3,2,1))

    println("Getting a small-scale displacement")
    CosmoMIA.filter_displacement!(alpt_storage.disp_field[1], args["disp_smooth"], box_size)
    CosmoMIA.filter_displacement!(alpt_storage.disp_field[2], args["disp_smooth"], box_size)
    CosmoMIA.filter_displacement!(alpt_storage.disp_field[3], args["disp_smooth"], box_size)

    

    patchy_field = UInt32.(round.(patchy_field))

    params = args["params"]

    Random.seed!(0)
    @time subgrid_cat = CosmoMIA.assign_particles_to_gals(alpt_storage, patchy_field, box_size, box_min, Normal(0, 0.6 * bin_size[1]); debug=false)
    
    @show reinterpret(SVector{3,eltype(subgrid_cat.pos)}, subgrid_cat.pos[:])[1]
    CosmoMIA.apply_displacement!(subgrid_cat.pos, alpt_storage.disp_field, args["step_size"], box_size, box_min)
    @show reinterpret(SVector{3,eltype(subgrid_cat.pos)}, subgrid_cat.pos[:])[1]


    alpt_storage = nothing
    
    collapsed_cat = CosmoMIA.copy_output(subgrid_cat)
    w = zeros(eltype(collapsed_cat.pos), size(collapsed_cat.pos, 2)) .+ 1

    function full_model(params, subgrid_cat, collapsed_cat; compute_power = true, compute_xi = true, do_collapse = true)

        collapsed_cat = CosmoMIA.copy_output(subgrid_cat)
        if do_collapse
            @time collapsed_cat = CosmoMIA.subgrid_collapse(collapsed_cat,
                                        params,
                                        box_size
                                        )
        end #if
        if compute_power
            
            collapsed_cat.pos = (collapsed_cat.pos .+ box_size) .% box_size
            k, power_real, _ = power_spectrum(collapsed_cat.pos, 
                                        w, 
                                        [0,2,4], 
                                        k_edges, 
                                        los, 
                                        grid_size, 
                                        box_size, 
                                        box_min, 
                                        2
                                        )
        else
            k = nothing
            power_real = nothing
        end #if
        collapsed_rsd = (CosmoMIA.apply_rsd!(zero(w), @view(collapsed_cat.pos[3,:]), @view(collapsed_cat.vel[3,:]), args["redshift"], cosmo))

        Threads.@threads for I = eachindex(collapsed_rsd)
            before = collapsed_rsd[I]
            if (collapsed_rsd[I] < 0) 
                collapsed_rsd[I] += box_size[3]
            elseif (collapsed_rsd[I] > box_size[3])
                collapsed_rsd[I] -= box_size[3]
                #println((collapsed_rsd[I] + box_size[3]) % box_size[3])
            end #if

            if (collapsed_rsd[I] < 0) | (collapsed_rsd[I] > box_size[3])
                error("Out of bounds coordinate after RSD for element $(collapsed_rsd[I]), before wrapping was $before")
            end #if
        end #for
        if compute_power
            k, power_rsd, _ = power_spectrum([@view(collapsed_cat.pos[1,:]), @view(collapsed_cat.pos[2,:]), collapsed_rsd],
                                        w, 
                                        [0,2,4], 
                                        k_edges, 
                                        los, 
                                        grid_size, 
                                        box_size, 
                                        box_min, 
                                        2
                                        )
        else
            power_rsd = nothing
        end #if
        if compute_xi
            s, xi_rsd = CosmoMIA.compute_2pcf(@view(collapsed_cat.pos[1,:]), @view(collapsed_cat.pos[2,:]), collapsed_rsd, box_size, s_edges)
        else
            s = nothing
            xi_rsd = nothing
        end #if
        k, power_real, power_rsd, s, xi_rsd                                     
    end #func

    k, power_real, power_rsd, s, xi_rsd = full_model(params, subgrid_cat, collapsed_cat, do_collapse = !args["skip_collapse"], compute_power = true, compute_xi = true)

    if args["plot"]
        p = plot!(p, k, k .^ 1.0 .* real.(@view(power_real[1, :])), xscale=:identity, ylabel=L"P(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Coll.")
        pp = plot(k, 100 .* (real.(@view(power_real[1, :])) ./ real.(@view(power_ref[1, :])) .- 1), label="Coll", ylabel=L"\frac{P}{P_{\mathrm{ref}}} -1 (\%)", ylim=(-5, 5))
        pp = plot!(pp, k, zero(k), ribbon=1, label=L"1\%", alpha=0.5)
    

        p0 = plot!(p0, k, k .^ 1.0 .* real.(@view(power_rsd[1, :])), xscale=:identity, ylabel=L"P(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Coll. l=0")
        p0 = plot!(p0, k, k .^ 1.0 .* real.(@view(power_rsd[2, :])), xscale=:identity, ylabel=L"P(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Coll. l=2")
        p0 = plot!(p0, k, k .^ 1.0 .* real.(@view(power_rsd[3, :])), xscale=:identity, ylabel=L"P(k)", xlabel=L"k~[h / \mathrm{Mpc}]", label="Coll. l=4")
        p0p = plot(k, 100 .* (real.(@view(power_rsd[1, :])) ./ real.(@view(power_ref_rsd[1, :])) .- 1), label="Coll", ylabel=L"\frac{P}{P_{\mathrm{ref}}} -1 (\%)", ylim=(-5, 5))
        p0p = plot!(p0p, k, 100 .* (real.(@view(power_rsd[2, :])) ./ real.(@view(power_ref_rsd[2, :])) .- 1), label="Coll", ylabel=L"\frac{P}{P_{\mathrm{ref}}} -1 (\%)", ylim=(-5, 5))
        p0p = plot!(p0p, k, zero(k), ribbon=1, label=L"1\%", alpha=0.5)

        px = plot!(px, s, s.^ 0 .* xi_rsd, label = "Coll,")
        pxp = plot(s, 100 .* (xi_rsd ./ xi_ref .- 1.))
        pxp = plot!(pxp, s, zero(s), ribbon=1, label=L"1\%", alpha=0.5, ylims = (-5, 5), xscale = :log10)
        final_plot = plot(p, p0, pp, p0p,px, pxp,
                      layout=@layout[
                                    a b
                                    c{0.1h} d{0.1h}
                                    _ e
                                    _ f{0.1h}
                                    ], 
                      dpi=300, 
                      size = (800, 800)
                      )
        savefig(final_plot, "plots/cosmoMIA_run.png")
    end #if
    if args["skip_collapse"]
        exit()
    end #if

    function loss(params; train_power = true, train_xi = true)
        if (params[3] < 0) || (params[4] < 0) || (params[5] < 0)
            return Inf
        end #if
        k, power_real, power_rsd, s, xi_rsd = full_model(params, subgrid_cat, collapsed_cat; compute_power = train_power, compute_xi = train_xi)
        if train_power
            mono_loss = mean(abs, @view(power_rsd[1,:]) .- @view(power_ref_rsd[1,:]))
            quad_loss = mean(abs, @view(power_rsd[2,:]) .- @view(power_ref_rsd[2,:]))
            hexa_loss = mean(abs, @view(power_rsd[3,:]) .- @view(power_ref_rsd[3,:]))
            real_loss = mean(abs, @view(power_real[1,:]) .- @view(power_ref[1,:]))
        else
            mono_loss = 0.
            quad_loss = 0.
            hexa_loss = 0.
            real_loss = 0.
        end #if
        if train_xi
            #xi_loss = mean(abs, s.^2 .* (xi_rsd .- xi_ref))
            xi_loss = mean(abs, (log10.(xi_rsd) .- log10.(xi_ref)))
        else
            xi_loss = 0.
        end #of
        total_loss = 3. * mono_loss + 1.5 * quad_loss + 0.5 * real_loss + 50. * xi_loss

        println("Total loss: $total_loss with params $params")
        
        total_loss
    end #func
    println("==> Starting optimization...")
    res = optimize((pars,) -> loss(pars; train_power = true, train_xi = true), params, NelderMead(), Optim.Options(show_trace = true))
    @show Optim.converged(res)
    params = Optim.minimizer(res)
    @show loss(Optim.minimizer(params))


end #func


function generate_from_params(args)
    cosmo = CosmoMIA.DESICosmology()
    box_size = @SVector [args["box_size"] for _ in 1:3]
    box_min = @SVector [0.0f0 for _ in 1:3]
    grid_size = Tuple(360 for _ in 1:3)
    k_edges = range(3f-3, 0.505, step = 5f-3)
    s_edges = Float32.(10 .^range(1f-1, log10(20), 51))
    los = [0,0,1]
    ref_size = 360
    bin_size = box_size ./ ref_size
    Random.seed!(42)

    println("Loading ALPT results")
    alpt_storage = CosmoMIA.ALPTResult{Float32}(box_size, ref_size; 
                                       pos_path = args["pos_path"], 
                                       vel_path = args["vel_path"], 
                                       disp_path = args["disp_path"], 
                                       delta_path = args["delta_path"],
                                       dweb_path = args["cweb_path"])
    println("Pre-boosting velocities.")
    CosmoMIA.renormalize_alpt_vel!(alpt_storage.vel_field, cosmo, args["redshift"])
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[1], args["vel_boost"], box_size)
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[2], args["vel_boost"], box_size)
    CosmoMIA.vel_kernel!(alpt_storage.vel_field[3], args["vel_boost"], box_size)
    
    

    patchy_field = permutedims(reshape(read!(args["number_counts"], Array{Float32}(undef, ref_size^3)), (ref_size, ref_size, ref_size)), (3,2,1))
    patchy_field = UInt32.(round.(patchy_field))

    params = args["params"]

    Random.seed!(0)
    @time subgrid_cat = CosmoMIA.assign_particles_to_gals(alpt_storage, patchy_field, box_size, box_min, Normal(0, 0.6 * bin_size[1]); debug=false)
    alpt_storage = nothing
    
    collapsed_cat = CosmoMIA.copy_output(subgrid_cat)
    w = zeros(eltype(collapsed_cat.pos), size(collapsed_cat.pos, 2)) .+ 1
    @time collapsed_cat = CosmoMIA.subgrid_collapse(collapsed_cat,
                                        params,
                                        box_size
                                        )

    #collapsed_rsd = (CosmoMIA.apply_rsd!(zero(w), @view(collapsed_cat.pos[3,:]), @view(collapsed_cat.vel[3,:]), args["redshift"], cosmo))


    npzwrite(args["output_path"], reduce(vcat, [collapsed_cat.pos, collapsed_cat.vel]))
    nothing

end #func

function main()
    args = parse_commandline()
    if args["do_fit"]
        fit_reference(args)
    else
        generate_from_params(args)
    end #if

end #func

main()