function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GswDrive,
    ind::Int,
    gsw::FT
) where {FT<:AbstractFloat}
    # update the conductances
    canopyi.g_sw[ind] = gsw
    gas_exchange!(photo_set, canopyi, envir, drive, ind)

    return nothing
end


function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GswDrive,
    ind::Int
) where {FT<:AbstractFloat}
    (; g_max, g_min) = canopyi

    # update the conductances
    canopyi.g_sw[ind] = clamp(canopyi.g_sw[ind], g_min, g_max)

    # update the conductances, 这里可以与其他模块合并
    canopyi.g_lw[ind] = 1 / (1 / canopyi.g_sw[ind] + 1 / canopyi.g_bw[ind])
    canopyi.g_sc[ind] = canopyi.g_sw[ind] / FT(1.6)
    canopyi.g_lc[ind] = 1 / (FT(1.6) / canopyi.g_sw[ind] + 1 / canopyi.g_m[ind] + 1 / canopyi.g_bc[ind])

    # update the photosynthetic rates
    leaf_photosynthesis!(photo_set, canopyi.ps, envir, GCO₂Mode(), canopyi.g_lc[ind])

    # be careful that this one might have memory allocation
    # make some changes on Photosynthesis.jl to avoid memory allocation
    # such as leaf_fluorescence!(photo_set, canopyi.ps);

    leaf_fluorescence!(photo_set.Flu, canopyi.ps)
    canopyi.Ac[ind] = canopyi.ps.Ac
    canopyi.Aj[ind] = canopyi.ps.Aj
    canopyi.Ap[ind] = canopyi.ps.Ap
    canopyi.Ag[ind] = canopyi.ps.Ag
    canopyi.An[ind] = canopyi.ps.An
    canopyi.φs[ind] = canopyi.ps.φs # Steady-state (light-adapted) yield

    # update the pressures
    canopyi.p_i[ind] = canopyi.ps.p_i # ps: Leaf photosynthesis system
    canopyi.p_s[ind] = canopyi.ps.p_s

    return nothing
end


function gas_exchange!(
    photo_set::AbstractPhotoModelParaSet{FT},
    canopyi::CanopyLayer{FT},
    envir::AirLayer{FT},
    drive::GswDrive
) where {FT<:AbstractFloat}
    # update the conductances for each "leaf"
    for i in eachindex(canopyi.g_sw)
        canopyi.ps.APAR = canopyi.APAR[i]
        leaf_ETR!(photo_set, canopyi.ps)
        gas_exchange!(photo_set, canopyi, envir, drive, i)
    end

    return nothing
end
