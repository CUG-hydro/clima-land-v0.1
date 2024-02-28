"""

    A_GROSS(node::SPACMono{FT}) where {FT}

Return the gross photosynthesis rate per layer, given
- `node` [`SPACMono`](@ref) type struct

"""
A_GROSS(node::SPACMono{FT}) where {FT} = [
    FT[iPS.Ag' * iPS.LAIx for iPS in node.plant_ps],
    FT[iPS.Ac' * iPS.LAIx for iPS in node.plant_ps],
    FT[iPS.Aj' * iPS.LAIx for iPS in node.plant_ps],
    FT[iPS.Ap' * iPS.LAIx for iPS in node.plant_ps]];


"""
    CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the canopy NPP of the SPAC per ground area
"""
function CNPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    cnpp::FT = 0
    for iPS in spac.plant_ps
        cnpp += numerical∫(iPS.An, iPS.LAIx) * iPS.LA
    end
    # cnpp = mapreduce(CNPP, +, spac.plant_ps)
    return cnpp / spac.ga
end;

"""
    GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}

Return the GPP of the SPAC per ground area
"""
function GPP(spac::SPACMono{FT}) where {FT<:AbstractFloat}
    gpp::FT = 0
    for iPS in spac.plant_ps
        gpp += numerical∫(iPS.Ag, iPS.LAIx) * iPS.LA
    end
    return gpp / spac.ga
end;



## For Canopy
function CNPP(iPS::CanopyLayer{FT})::FT where {FT<:AbstractFloat}
    numerical∫(iPS.An, iPS.LAIx) * iPS.LA # [μmol m⁻² s⁻¹] * [m²]
end

function GPP(iPS::CanopyLayer{FT})::FT where {FT<:AbstractFloat}
    numerical∫(iPS.Ag, iPS.LAIx) * iPS.LA # [μmol m⁻² s⁻¹] * [m²]
end
