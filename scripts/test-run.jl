using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using PkgUtility
using Test

FT = Float64
AirLayer{FT}()
Leaf{FT}()

C3CLM(FT)
C4CLM(FT)

envir = AirLayer{FT}();
leaf_3 = Leaf{FT}();
leaf_4 = Leaf{FT}();
mod_3 = C3CLM(FT);
mod_4 = C4CLM(FT);
rand_T = rand(FT) + 298;

can_3 = CanopyLayer{FT}(n_leaf=2);
can_4 = CanopyLayer{FT}(n_leaf=2);
esm_1 = ESMBallBerry{FT}();
esm_2 = ESMGentine{FT}();
esm_3 = ESMLeuning{FT}();
esm_4 = ESMMedlyn{FT}();
osm_1 = OSMEller{FT}();
osm_2 = OSMSperry{FT}();
osm_3 = OSMWang{FT}();
osm_4 = OSMWAP{FT}();
osm_5 = OSMWAPMod{FT}();
hs = LeafHydraulics{FT}();

# Update leaf physiological parameters if temperature or pressure changes in the daytime
update_leaf_TP!(mod_3, can_3, hs, envir);

# # Update leaf maximal A and K for Sperry model, can skip
# update_leaf_AK!(mod_3, can_3, hs, envir);

gas_exchange!(mod_3, can_3, envir, GlcDrive(), 1, FT(0.1));
gas_exchange!(mod_3, can_3, envir, GswDrive(), 1, FT(0.05));
# @test PkgUtility.NaN_test(can_3);
# @test PkgUtility.NaN_test(can_4);

can_3.g_sw[2] = 0;
gsw_control!(mod_3, can_3, envir, 2);
# @test PkgUtility.NaN_test(can_3);
# @test PkgUtility.NaN_test(can_4);
