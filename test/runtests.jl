using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using PkgUtility
using Test

using Land.SoilPlantAirContinuum: CNPP
FT = Float32



ENV["JULIA_LOG_LEVEL"] = "WARN"

@testset verbose = true "CliMA Land v0.1" begin
    include("modules/CanopyLayers.jl")
    include("modules/Photosynthesis.jl")
    include("modules/PlantHydraulics.jl")
    include("modules/StomataModels.jl")
    include("modules/SPAC.jl")
end;

@testset verbose = true "CliMA Land Features" begin
    include("features/clm5_mode.jl")
end;


@testset "CNPP, GPP, T_VEG" begin
    canopy = CanopyLayer{FT}()
    envir = AirLayer{FT}()

    CNPP(canopy)
    GPP(canopy)
    T_VEG(canopy, envir)
end
