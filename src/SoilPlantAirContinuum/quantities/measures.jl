# 1. 冠层分层
# 2. 每层很多叶子

"""
    T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the transpiration of the SPAC per ground area
"""
function T_VEG(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    t_veg::FT = 0;

    for i_can in 1:spac.n_canopy
        iEN = spac.envirs[i_can];
        iPS = spac.plant_ps[i_can];
        t_veg += numerical∫(iPS.g_lw, iPS.LAIx) * (iPS.p_sat - iEN.p_H₂O) / iEN.p_atm * iPS.LA;
    end;

    return t_veg / spac.ga # [mol s⁻¹] / [m⁻²] = [mol m⁻² s⁻¹]
end


function T_VEG(plant::CanopyLayer, envir::AirLayer{FT}) where {FT<:AbstractFloat}
    (; g_lw, LAIx, LA, p_sat) = plant
    (; p_H₂O, p_atm) = envir
    # LA: Total leaf area [m²]
    # gs: water vapor conductance [mol s⁻¹] at the canopy level
    gw = numerical∫(g_lw, LAIx) * LA # [mol m⁻² s⁻¹] * [m²] = [mol s⁻¹]
    trans = gw * (p_sat - p_H₂O) / p_atm # T_VEG, [mol s⁻¹]
    # trans * M_H₂O() # kg s⁻¹, use unit object
    # / spac.ga # [mol s⁻¹] * [kg mol⁻¹] / [m²] = [kg m⁻² s⁻¹]
end


"""
    PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the cumulative PPAR of the SPAC per ground area
"""
function PPAR(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    ppar::FT = 0;

    for iPS in spac.plant_ps
        ppar += numerical∫(iPS.APAR, iPS.LAIx) * FT(1e-6) * spac.canopy_rt.LAI / spac.canopy_rt.nLayer;
    end;

    return ppar
end
