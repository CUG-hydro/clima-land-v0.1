function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GlcDrive,
    ind::Int,
    glc::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_lc[ind] = glc
    gas_exchange!(photo_set, canopyi, envir, drive, ind)

    return nothing
end



function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GlcDrive,
    ind::Int
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sc[ind] = 1 / (1 / canopyi.g_lc[ind] - 1 / canopyi.g_m[ind] - 1 / canopyi.g_bc[ind])
    canopyi.g_sw[ind] = canopyi.g_sc[ind] * FT(1.6)
    canopyi.g_lw[ind] = 1 / (1 / canopyi.g_sw[ind] + 1 / canopyi.g_bw[ind])

    # update the photosynthetic rates
    if canopyi.g_lc[ind] != canopyi.ps.g_lc
        leaf_photosynthesis!(photo_set, canopyi.ps, envir, GCO₂Mode(), canopyi.g_lc[ind])

        # be careful that this one might have memory allocation
        # make some changes on Photosynthesis.jl to avoid memory allocation
        # such as leaf_fluorescence!(photo_set, canopyi.ps);
        leaf_fluorescence!(photo_set.Flu, canopyi.ps)
    end
    canopyi.Ac[ind] = canopyi.ps.Ac
    canopyi.Aj[ind] = canopyi.ps.Aj
    canopyi.Ap[ind] = canopyi.ps.Ap
    canopyi.Ag[ind] = canopyi.ps.Ag
    canopyi.An[ind] = canopyi.ps.An
    canopyi.φs[ind] = canopyi.ps.φs

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_i
    canopyi.p_s[ind] = canopyi.ps.p_s

    return nothing
end




function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GlcDrive
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_lc)
        canopyi.ps.APAR = canopyi.APAR[i]
        leaf_ETR!(photo_set, canopyi.ps)
        gas_exchange!(photo_set, canopyi, envir, drive, i, canopyi.g_lc[i])
    end

    return nothing
end
