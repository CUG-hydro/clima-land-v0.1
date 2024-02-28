"""

    run_time_step!(spac::SPACMono{FT}, dfr::DataFrame) where {FT<:AbstractFloat}

Run CliMA Land in a time step, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `ind` Time index

"""
function run_time_step!(spac::SPACMono{FT}, dfr::DataFrameRow, beta::BetaGLinearPsoil{FT};
    use_sif=true) where {FT<:AbstractFloat}

    # read the data out of dataframe row to reduce memory allocation
    _df_dif::FT = dfr.RAD_DIF
    _df_dir::FT = dfr.RAD_DIR

    # compute beta factor (based on Psoil, so canopy does not matter)
    _βm = spac_beta_max(spac, beta)

    # calculate leaf level flux per canopy layer
    for _i_can in 1:spac.n_canopy
        _iPS = spac.plant_ps[_i_can] # CanopyLayer
        _iEN = spac.envirs[_i_can]   # env

        # set gsw to 0 or iterate for 30 times to find steady state solution
        if _df_dir + _df_dif < 10
            # night
            _iPS.APAR .= 0
            _iPS.g_sw .= 0 # Stomatal conductance to water H₂O
            gsw_control!(spac.photo_set, _iPS, _iEN)
        else
            for _ in 1:30
                gas_exchange!(spac.photo_set, _iPS, _iEN, GswDrive())
                update_gsw!(spac, spac.stomata_model, _i_can, FT(120); β=_βm) # 求解gsw, 更新叶片导度
                # prognostic_gsw!(spac.plant_ps[_i_can], spac.envirs[_i_can], spac.stomata_model, βm, FT(120))
                gsw_control!(spac.photo_set, _iPS, _iEN)
            end
        end
    end

    # calculate the SIF if there is sunlight
    if _df_dir + _df_dif >= 10 && use_sif
        # 不需要则不计算
        update_sif!(spac)
        dfr.SIF683 = SIF_WL(spac.can_rad, spac.wl_set, FT(682.5))
        dfr.SIF740 = SIF_740(spac.can_rad, spac.wl_set)
        dfr.SIF757 = SIF_WL(spac.can_rad, spac.wl_set, FT(758.7))
        dfr.SIF771 = SIF_WL(spac.can_rad, spac.wl_set, FT(770.0))
        dfr.NDVI = NDVI(spac.can_rad, spac.wl_set)
        dfr.EVI = EVI(spac.can_rad, spac.wl_set)
        dfr.NIRv = NIRv(spac.can_rad, spac.wl_set)
    end

    # save the total flux into the DataFrame
    dfr.F_H2O = T_VEG(spac) # transpiration
    dfr.F_CO2 = CNPP(spac) # net photosynthesis
    dfr.F_GPP = GPP(spac)

    return nothing
end


"""

    run_model!(spac::SPACMono{FT}, df::DataFrame, nc_out::String) where {FT<:AbstractFloat}

Run CliMA Land at a site for the enture year, given
- `spac` Soil plant air continuum struct
- `df` Weather driver dataframe
- `nc_out` File path to save the model output

"""
function run_model!(spac::SPACMono{FT}, df::DataFrame; use_sif=true) where {FT<:AbstractFloat}
    _in_rad_bak = deepcopy(spac.in_rad)
    _in_dir = _in_rad_bak.E_direct' * spac.wl_set.dWL / 1000
    _in_dif = _in_rad_bak.E_diffuse' * spac.wl_set.dWL / 1000
    _deepcopies = [_in_rad_bak, _in_dir, _in_dif]
    _beta_g = BetaGLinearPsoil{FT}()

    # set up memory
    _spac_mem = SPACMemory{FT}()

    # iterate through the time steps
    # add a progress bar to show the progress
    @showprogress for _dfr in eachrow(df)
        prescribe_parameters!(spac, _dfr, _spac_mem, _deepcopies)
        run_time_step!(spac, _dfr, _beta_g; use_sif)
    end

    # # save simulation results to hard drive
    # save_nc!(nc_out, df[:, DF_VARIABLES])
    df[:, DF_VARIABLES]
    # return nothing
end
