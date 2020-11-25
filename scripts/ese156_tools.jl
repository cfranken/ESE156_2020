using Pkg.Artifacts
using NCDatasets
using ProgressMeter
using Plots
using ImageFiltering
using Distributions
using Interpolations
using Polynomials

using RadiativeTransfer.CrossSection
using Parameters
using CSV
using DelimitedFiles
using ForwardDiff, DiffResults

abstract type InstrumentOperator end

"Struct for an atmospheric profile"
struct AtmosphericProfile{FT}
    lat::FT
    lon::FT
    psurf::FT
    T::Array{FT,1}
    q::Array{FT,1}
    p::Array{FT,1}
    p_levels::Array{FT,1}
    vmr_h2o::Array{FT,1}
    vcd_dry::Array{FT,1}
    vcd_h2o::Array{FT,1}
end;

"Struct for Kernel Instrument Function"
@with_kw struct KernelInstrument{FT} <: InstrumentOperator 
    kernel::FT
    ν_out::Array
end;

"Read atmospheric profile (just works for our file, can be generalized"
function read_atmos_profile(file::String, lat::Real, lon::Real, timeIndex; g₀=9.8196)
    @assert 1 <= timeIndex <= 4

    ds = Dataset(file)

    # See how easy it is to actually extract data? Note the [:] in the end reads in ALL the data in one step
    lat_   = ds["YDim"][:]
    lon_   = ds["XDim"][:]
    
    FT = eltype(lat_)
    lat = FT(lat)
    lon = FT(lon)
    
    # Find index (nearest neighbor, one could envision interpolation in space and time!):
    iLat = argmin(abs.(lat_ .- lat))
    iLon = argmin(abs.(lon_ .- lon))
    @show ds["T"]
    # Temperature profile
    T    = convert(Array{FT,1}, ds["T"][ iLon,iLat, :, timeIndex])
    # specific humidity profile
    q    = convert(Array{FT,1}, ds["QV"][iLon,iLat,  :, timeIndex])
    
    # Surafce pressure
    psurf = convert(FT, ds["PS"][iLon, iLat, timeIndex])
    
    # AK and BK global attributes (important to calculate pressure half-levels)
    ak = ds.attrib["HDF_GLOBAL.ak"][:]
    bk = ds.attrib["HDF_GLOBAL.bk"][:]

    p_half = (ak + bk * psurf)
    p_full = (p_half[2:end] + p_half[1:end - 1]) / 2
    close(ds)
    
    # Avogradro's number:
    Na = 6.0221415e23;
    # Dry and wet mass
    dryMass = 28.9647e-3  / Na  # in kg/molec, weighted average for N2 and O2
    wetMass = 18.01528e-3 / Na  # just H2O
    ratio = dryMass / wetMass 
    n_layers = length(T)
    # also get a VMR vector of H2O (volumetric!)
    vmr_h2o = zeros(FT, n_layers, )
    vcd_dry = zeros(FT, n_layers, )
    vcd_h2o = zeros(FT, n_layers, )

    # Now actually compute the layer VCDs
    for i = 1:n_layers 
        Δp = p_half[i + 1] - p_half[i]
        vmr_h2o[i] = q[i] * ratio
        vmr_dry = 1 - vmr_h2o[i]
        M  = vmr_dry * dryMass + vmr_h2o[i] * wetMass
        vcd_dry[i] = vmr_dry * Δp / (M * g₀ * 100.0^2)   # includes m2->cm2
        vcd_h2o[i] = vmr_h2o[i] * Δp / (M * g₀ * 100^2)
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_half, vmr_h2o, vcd_dry, vcd_h2o)
end;

"Computes cross section matrix for arbitrary number of absorbers"
function compute_profile_crossSections(profile::AtmosphericProfile, hitranModels::Array{HitranModel}, ν::AbstractRange{<:Real})
    nGases   = length(hitranModels)
    nProfile = length(profile.p)
    FT = eltype(profile.T)
    n_layers = length(profile.T)

    σ_matrix = zeros(FT, (length(ν), n_layers, nGases))
    @showprogress for i = 1:n_layers
        p_ = profile.p[i] / 100 # in hPa
        T_ = profile.T[i]
        for j = 1:nGases
            σ_matrix[:,i,j] = absorption_cross_section(hitranModels[j], ν, p_, T_);
        end
    end
    return σ_matrix
end;

"Creates a Gaussian Kernel"
function gaussian_kernel(FWHM, res; fac=5)
    co = 2.355
    width = FWHM / res / co
    extent = ceil(fac * width)
    d = Normal(0, width)
    kernel = centered(pdf.(d, -extent:extent))
    return kernel ./ sum(kernel)
end;

"Rescales x-axis to go from [-1,1]"
function rescale_x(a)
    a = a .- mean(a);
    a = a / (a[end] - a[1]) * 2;
end;

"Convolve and resample"
function conv_spectra(m::KernelInstrument, ν, spectrum)
    s = imfilter(spectrum, m.kernel)
    interp_cubic = CubicSplineInterpolation(ν, s)
    return interp_cubic(m.ν_out)
end;
    

"Reduce profile dimensions"
function reduce_profile(n::Int, profile::AtmosphericProfile, σ_matrix)
    @assert n < length(profile.T)
    @unpack lat, lon, psurf = profile
    # New rough half levels (boundary points)
    a = range(0, maximum(profile.p), length=n + 1)
    dims = size(σ_matrix)
    FT = eltype(σ_matrix)
    σ_matrix_lr = zeros(FT, dims[1], n, dims[3])
    T = zeros(FT, n);
    q = zeros(FT, n);
    p_full = zeros(FT, n);
    p_levels = zeros(FT, n + 1);
    vmr_h2o  = zeros(FT, n);
    vcd_dry  = zeros(FT, n);
    vcd_h2o  = zeros(FT, n);

    for i = 1:n
        ind = findall(a[i] .< profile.p .<= a[i + 1]);
        σ_matrix_lr[:,i,:] = mean(σ_matrix[:,ind,:], dims=2);
        p_levels[i] = profile.p_levels[ind[1]]
        p_levels[i + 1] = profile.p_levels[ind[end]]
        p_full[i] = mean(profile.p_levels[ind])
        T[i] = mean(profile.T[ind])
        q[i] = mean(profile.q[ind])
        vmr_h2o[i] = mean(profile.vmr_h2o[ind])
        vcd_dry[i] = sum(profile.vcd_dry[ind])
        vcd_h2o[i] = sum(profile.vcd_h2o[ind])
    end

    return AtmosphericProfile(lat, lon, psurf, T, q, p_full, p_levels, vmr_h2o, vcd_dry, vcd_h2o), σ_matrix_lr
end;

function artifact_ESE_folder(uuid::Union{String,Nothing}=nothing)
    name = "ESE156_2020"
    artifacts_toml = find_artifacts_toml(pwd())
    artifact_dict = Artifacts.load_artifacts_toml(artifacts_toml)
    Artifacts.do_artifact_str(name, artifact_dict, artifacts_toml, uuid) 
end


    
    
    




    
