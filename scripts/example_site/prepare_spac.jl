"""

    prepare_spac(dict::Dict; FT = Float64)

Create a SPAC, given
- `dict` Dictionary of GriddingMachine data in a grid
- `FT` Floating number type

"""
function prepare_spac(dict::Dict; FT=Float64)
    # read general information from dict
    _lat = dict["latitude"]
    _lon = dict["longitude"]
    _sm = ESMMedlyn{FT}()

    # use JULES soil depth 0.00 -- 0.10 -- 0.35 -- 1.00 -- 3.00 m, and assume 2 m deep root (z_root = -2) for all the sites
    soil_bounds = FT[0, -0.1, -0.35, -1, -3]
    z_canopy = max(FT(0.1), dict["canopy_height"])
    Δz = z_canopy / 20
    air_bounds = collect(0:Δz:z_canopy+2*Δz)
    plant_hs = create_tree(FT(-2), z_canopy / 2, z_canopy, soil_bounds, air_bounds)

    # create a SPACMono struct, redefine the wavelength limits for PAR if ePAR is true
    _node = SPACMono{FT}(;soil_bounds, air_bounds, z_canopy, z_root=-2, plant_hs, 
        latitude=_lat, longitude=_lon, stomata_model=_sm)

    for _iPS in _node.plant_ps
        _iPS.g_min = eps(FT)
        _iPS.g_min25 = eps(FT)
        _iPS.g_max = 0.8
        _iPS.g_max25 = 0.8
    end

    # update soil type information per layer
    for _i in eachindex(_node.plant_hs.roots)
        α = dict["soil_vg_α"][_i]
        n = dict["soil_vg_n"][_i]
        Θr = dict["soil_vg_Θr"][_i]
        Θs = dict["soil_vg_Θs"][_i]
        _node.plant_hs.roots[_i].sh = VanGenuchten{FT}(;stype="JULES", α, n, Θs, Θr)
    end

    # update leaf mass per area (from m² kg⁻¹ to g cm⁻²)
    _lma = dict["leaf_mass_per_area"]
    for leaf in _node.leaves_rt
        leaf.Cm = _lma
    end

    # set up empirical model
    if typeof(_sm) <: ESMMedlyn
        _node.photo_set = C3CLM(FT)
        _node.stomata_model.g1 = dict["g1_medlyn_c3"]
        _node.stomata_model.g0 = 1e-3
    else
        @warn "Stomatal model parameters are not initialized for $(typeof(_sm))"
    end

    # update soil color class from CLM dataset
    _node.soil_opt.color = dict["soil_color"]

    # update the Vcmax, Jmax, and Vpmax
    update_VJRWW!(_node, nanmean(dict["vcmax"]))

    # initialize the canopy RT model
    initialize_spac_canopy!(_node)

    return _node
end
