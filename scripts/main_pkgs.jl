using DataFrames: DataFrame, DataFrameRow
using Dates: isleapyear
using JLD2: load

using Revise
using ProgressMeter
using UnPack

using NetcdfIO: read_nc, save_nc!
using PkgUtility: month_days, nanmean

using Land.CanopyLayers: EVI, FourBandsFittingHybrid, NDVI, NIRv, SIF_WL, SIF_740, fit_soil_mat!
using Land.Photosynthesis: C3CLM, use_clm_td!
using Land.PlantHydraulics: VanGenuchten, create_tree
using Land.SoilPlantAirContinuum: CNPP, GPP, PPAR, SPACMono, T_VEG, initialize_spac_canopy!, 
    prescribe_air!, prescribe_swc!, prescribe_t_leaf!, spac_beta_max, update_Cab!, update_LAI!, update_VJRWW!,
    update_par!, update_sif!, zenith_angle
using Land.StomataModels: BetaGLinearPsoil, ESMMedlyn, GswDrive, gas_exchange!, gsw_control!, prognostic_gsw!

using Land
using Land.CanopyLayers
using Land.Photosynthesis
using Land.PlantHydraulics
using Land.SoilPlantAirContinuum
using Land.StomataModels
using PkgUtility
using Test

# DF_VARIABLES = ["F_H2O", "F_CO2", "F_GPP", "SIF683", "SIF740", "SIF757", "SIF771", "NDVI", "EVI", "NIRv"];
DF_VARIABLES = ["F_H2O", "F_CO2", "F_GPP"];

# Structure that store memory information
Base.@kwdef mutable struct SPACMemory{FT<:AbstractFloat}
    chl::FT = -9999
    lai::FT = -9999
    vcm::FT = -9999
end

# include("main_input.jl")
includet("./example_site/run_time_step!.jl")
includet("./example_site/prescribe_parameters!.jl")
includet("./example_site/prepare_wd.jl")
includet("./example_site/prepare_spac.jl")
