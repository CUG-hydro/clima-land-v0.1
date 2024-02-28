FT = Float32
c3_set = C3CLM(FT)
# c4_set = C4CLM(FT)
leaf = Leaf{FT}()
# leaf_4 = Leaf{FT}()
td_q10 = Photosynthesis.Q10TD{FT}(1, 273.15, 1.7)
envir = AirLayer{FT}()
fluo_set = c3_set.Flu
T = rand(FT) + 298
glc = FT(0.1)
p_i = rand(FT) + 20

# temperature corrections
Photosynthesis.photo_TD_from_set(td_q10, T)
leaf_temperature_dependence!(c3_set, leaf, envir, T)

# rubisco limited rates
Photosynthesis.rubisco_limited_rate!(c3_set, leaf)
Photosynthesis.rubisco_limited_rate!(c3_set, leaf, envir)
# @test PkgUtility.NaN_test(leaf)
# @test PkgUtility.NaN_test(leaf)

# light limited rates
Photosynthesis.leaf_ETR!(c3_set, leaf)
Photosynthesis.light_limited_rate!(c3_set, leaf)
# @test PkgUtility.NaN_test(leaf)
Photosynthesis.light_limited_rate!(c3_set, leaf, envir)
# @test PkgUtility.NaN_test(leaf)

# product limited rates
Photosynthesis.product_limited_rate!(c3_set, leaf)
@test PkgUtility.NaN_test(leaf)

# fluorescence
leaf_photosynthesis!(c3_set, leaf, envir, PCO₂Mode(), FT(2))
leaf_fluorescence!(fluo_set, leaf)

leaf_photosynthesis!(c3_set, leaf, envir, GCO₂Mode())
leaf_fluorescence!(fluo_set, leaf)
# @test PkgUtility.NaN_test(leaf)

# leaf photo from glc
leaf_photosynthesis!(c3_set, leaf, envir, GCO₂Mode(), glc)
# @test PkgUtility.NaN_test(leaf)

# leaf photo from p_i
leaf_photosynthesis!(c3_set, leaf, envir, PCO₂Mode(), p_i)
# @test PkgUtility.NaN_test(leaf)
