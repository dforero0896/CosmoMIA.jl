
function ∂χ_flat_LCDM(z, H0, Ω_Λ, Ω_m, speed_of_light = 299792.458)
    H = H0 * √(Ω_Λ + Ω_m * (1 + z)^3)
    speed_of_light / H
end #func
function χ_flat_LCDM(z, H0, Ω_Λ, Ω_m, speed_of_light = 299792.458)
    quadgk(z -> ∂χ_flat_LCDM(z, H0, Ω_Λ, Ω_m), 0., z, rtol=1e-8)[1]
end #func


const eV = 1.782661907e-36  # [kg] - eV/c^2
const mpc = 3.085677581491367e22  # [m]
const speed_of_light_km_s = 299792.458
const G = 6.67408e-11  # [m^3/kg/s^2]
const aRad = 4.0*5.670400e-8 # radiation constant
const σT = 6.6524587158e-29  # [m**2]
const mp = 938.2720813  # [MeV/c**2]
const msun = 1.98848e30  # [kg]
const ħ = 6.582119514e-16  # eV s
const evc² = 1.782661907e-36  # [kg] - eV/c^2
const kb = 8.6173303e-5 
const cgs_Mpc = 1.
const cgs_sec = 1.
const cgs_km = 0.3240779290e-19 * cgs_Mpc
const cgs_clight = (0.9715611892e-14 * cgs_Mpc/cgs_sec)
###
#define cgs_Mpc static_cast<real_prec>(1.)
#define cgs_sec static_cast<real_prec>(1.)
#define cgs_km static_cast<real_prec>((0.3240779290e-19 * cgs_Mpc))
#define cgs_clight static_cast<real_prec>((0.9715611892e-14 * cgs_Mpc/cgs_sec))
###


@with_kw mutable struct Cosmology

    ω_b = 0.0225f0
    ω_c = 0.12f0
    h = 0.67f0
    Neff = 3.044f0
    Ω_k₀ = 0.f0
    T_cmb = 2.725f0
    w0 = -1f0
    wa = 0f0
    ns = 0.96f0
    T = typeof(ω_b)
    h² = T(h^2)
    H₀ = h * 100
    
    # Photons
    T_γ₀ = T_cmb #* Kelvin
    ρ_γ₀ = T(3 * 100^2 / (8π * G) * ħ^3 * (1. / (1e-3 * mpc))^2 * (1. / evc²) * (speed_of_light_km_s * 1e3)^3)
    Ω_γ₀ = T((π^2 / 15) * (2.725 * kb)^4 / ρ_γ₀ * (T_cmb / 2.725)^4 / h²)
    ω_γ = Ω_γ₀ * h²
    # Neutrinos
    Ω_ν₀ = T(Neff * (7 / 8 * (4 / 11)^(4 / 3)) * Ω_γ₀)
    ω_ν = T(Ω_ν₀ * h²)
    T_ν₀ = T_γ₀ * (4 / 11)^(1 / 3)
    # baryons
    Ω_b₀ = ω_b / h²
    # CDM
    Ω_c₀ = ω_c / h²
    # Curvature
    ω_k = Ω_k₀ * h²
    # Dark energy
    Ω_Λ₀ = 1 - Ω_k₀ - Ω_b₀ - Ω_c₀ - Ω_ν₀ - Ω_γ₀
    ω_Λ = Ω_Λ₀ * h²

    Ω_m₀ = Ω_b₀ + Ω_c₀
    Ω_r₀ = Ω_γ₀ + Ω_ν₀

    cache = nothing
    z_tab_min = 0
    z_tab_max = 3
    z_tab_num = 100000

end #struct

function DESICosmology(;kwargs...)
    Cosmology(;ω_b =0.02237f0, ω_c = 0.1200f0, h = 0.6736f0, ns = 0.9649f0, Neff = 2.0328f0 + 1f0, w0 = -1f0, wa = 0f0, kwargs...)
end #funct

Ω_ν(c::Cosmology, z) = c.Ω_ν₀ * (1 + z)^4
Ω_γ(c::Cosmology, z) = c.Ω_γ₀ * (1 + z)^4
Ω_b(c::Cosmology, z) = c.Ω_b₀ * (1 + z)^3
Ω_c(c::Cosmology, z) = c.Ω_c₀ * (1 + z)^3
Ω_k(c::Cosmology, z) = c.Ω_k₀ * (1 + z)^2
w(c::Cosmology, z) = c.w0 + (1 - (1+z)^-1) * c.wa
f_DE(c::Cosmology, a) = (-3 * (1 + c.w0) + 3 * c.wa * ((a - 1) / log(a) - 1))
Ω_Λ(c::Cosmology, z) = c.Ω_Λ₀ * (1 + z)^(-f_DE(c, (1+z)^-1))
E(c::Cosmology, z) = sqrt(Ω_ν(c, z) + Ω_γ(c, z) + Ω_b(c,z) + Ω_c(c,z) + Ω_k(c,z) +  Ω_Λ(c, z))
H(c::Cosmology, z) = c.h * 100 * E(c, z)
χ(c::Cosmology, z1, z2) = quadgk(z -> 1 / H(c, z), z1, z2, rtol=1e-8)[1] # Mpc s / km
χ(c::Cosmology, z) = χ(c, 0, z) # Mpc s / km
comoving_distance(c::Cosmology, z) = speed_of_light_km_s * χ(c, z)  # Mpc
function comoving_distance_interp(c::Cosmology)
    if c.cache === nothing
        c.cache = Dict()
        c.cache[:z] = range(c.z_tab_min, c.z_tab_max, length=c.z_tab_num)
        c.cache[:r] = comoving_distance.(Ref(c), c.cache[:z])
    end #if
    interpolate((c.cache[:z],), c.cache[:r], Gridded(Linear()))
end #func
function redshift_interp(c::Cosmology)
    if c.cache === nothing
        c.cache = Dict()
        c.cache[:z] = range(c.z_tab_min, c.z_tab_max, length=c.z_tab_num)
        c.cache[:r] = comoving_distance.(Ref(c), c.cache[:z])
    end #if
    interpolate((c.cache[:r],), c.cache[:z], Gridded(Linear()))
end #func

function growth_rate_approx(c::Cosmology, z; order = 1)
    # Eq. 11a/b https://ui.adsabs.harvard.edu/abs/1991MNRAS.251..128L/abstract
    Ω = c.Ω_m₀ * E(c, z)^-2 * (1 + z)^3
    if order == 1
        return Ω^(5. / 9)
    elseif order == 2
        return 2. * Ω^(6. /11)
    else
        error("Growth rate implemented to second order")
    end #if

end #func
#function growth_rate_approx(c::Cosmology, z)
#    # Eq. 11a/b https://ui.adsabs.harvard.edu/abs/1991MNRAS.251..128L/abstract
#    Ω = c.Ω_m₀ * E(c, z)^-2 * (1 + z)^3
#    Ω^(5. / 9)
#end #func
function growth_factor_approx(cosmo::Cosmology, z)
    a = 1. / 3
    b = 1.
    c = 11. / 6
    d = (1+z)^(-3) * (1. - 1. / cosmo.Ω_m₀)
    HypergeometricFunctions._₂F₁(a, b, c, d) / (1 + z)
end #func