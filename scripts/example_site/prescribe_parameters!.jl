"""

    prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrame, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}

Prescibe parameters for the SPAC, given
- `spac` Soil plant air continuum struct
- `dfr` Weather driver dataframe row
- `mem` Memory cache struct
- `deepcopies` Deepcopies of radiation used to scale direct and diffuse radiation

"""
function prescribe_parameters!(spac::SPACMono{FT}, dfr::DataFrameRow, mem::SPACMemory{FT}, deepcopies::Vector) where {FT<:AbstractFloat}
    # read the data out of dataframe row to reduce memory allocation
    _df_atm::FT = dfr.P_ATM
    _df_chl::FT = dfr.Chlorophyll
    _df_cli::FT = dfr.CI
    _df_co2::FT = dfr.CO2
    _df_dif::FT = dfr.RAD_DIF
    _df_dir::FT = dfr.RAD_DIR
    _df_doy::FT = dfr.FDOY
    _df_lai::FT = dfr.LAI
    _df_sw1::FT = dfr.SWC_1
    _df_sw2::FT = dfr.SWC_2
    _df_sw3::FT = dfr.SWC_3
    _df_sw4::FT = dfr.SWC_4
    _df_tar::FT = dfr.T_AIR
    _df_tlf::FT = dfr.T_LEAF
    _df_tmn::FT = dfr.T_MEAN
    _df_vcm::FT = dfr.Vcmax
    _df_vpd::FT = dfr.VPD
    _df_wnd::FT = dfr.WIND

    # adjust optimum t based on 10 day moving average skin temperature
    use_clm_td!(spac.photo_set, _df_tmn)

    # if total LAI, Vcmax, or Chl changes, update them (add vertical Vcmax profile as well)
    _trigger_lai::Bool = !isnan(_df_lai) && (_df_lai != mem.lai)
    _trigger_vcm::Bool = !isnan(_df_vcm) && (_df_vcm != mem.vcm)
    _trigger_chl::Bool = !isnan(_df_chl) && (_df_chl != mem.chl)
    if _trigger_lai
        update_LAI!(spac, _df_lai)
        mem.lai = _df_lai
    end

    if _trigger_lai || _trigger_vcm
        update_VJRWW!(spac, _df_vcm; expo=FT(0.3))
        mem.vcm = _df_vcm
    end

    if _trigger_chl
        update_Cab!(spac, _df_chl; cab_2_car=FT(1 / 7))
        mem.chl = _df_chl
    end

    # update clumping index
    spac.canopy_rt.Ω = _df_cli
    spac.canopy_rt.clump_a = _df_cli

    # sync the environmental conditions per layer
    prescribe_air!(spac, _df_co2, _df_atm, _df_tar, _df_vpd, _df_wnd)
    prescribe_t_leaf!(spac, max(_df_tar, _df_tlf))

    # run the chunks below only when total radiation is higher than 10
    if _df_dir + _df_dif < 10
        return nothing
    end

    # update soil water matrices per layer
    prescribe_swc!(spac, _df_sw1, _df_sw2, _df_sw3, _df_sw4)

    # update soil albedo using FourBandsFittingHybrid
    _method = FourBandsFittingHybrid()
    fit_soil_mat!(spac.soil_opt, spac.wl_set, spac.swc[1], _method)

    # update PAR related information
    spac.in_rad.E_direct .= deepcopies[1].E_direct .* _df_dir ./ deepcopies[2]
    spac.in_rad.E_diffuse .= deepcopies[1].E_diffuse .* _df_dif ./ deepcopies[3]
    spac.angles.sza = min(88, zenith_angle(spac.latitude, _df_doy))
    update_par!(spac)

    return nothing
end
