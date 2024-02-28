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

# 这里为何没有考虑冠层辐射传输？
canopy = CanopyLayer{FT}()
envir = AirLayer{FT}()

CNPP(canopy)
T_VEG(envir, canopy)
