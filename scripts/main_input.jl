"""

    prepare_wd(dict::Dict, wd_file::String)

Prepare weather driver dataframe to feed CliMA Land, given
- `dict` Dictionary that store grid information
- `wd_file` Weather driver file

"""
function prepare_wd(dict::Dict, wd_file::String)
  _df_in = read_nc(wd_file)

  # compute T_MEAN based on the weather driver
  _df_in[!, "CO2"] .= 0.0
  _df_in[!, "Chlorophyll"] .= 0.0
  _df_in[!, "LAI"] .= 0.0
  _df_in[!, "Vcmax"] .= 0.0
  _df_in[!, "T_MEAN"] .= 0.0
  for _i in eachindex(_df_in.T_MEAN)
    if _i < 240
      _df_in[_i, "T_MEAN"] = nanmean(max.(_df_in.T_AIR[1:_i], _df_in.T_LEAF[1:_i]))
    else
      _df_in[_i, "T_MEAN"] = nanmean(max.(_df_in.T_AIR[_i-239:_i], _df_in.T_LEAF[_i-239:_i]))
    end
  end

  #
  # extropolate the time series based on input variable dimensions, dimensions must be within supported settings
  #     1. extropolate the data to 1D resolution
  #     2. extropolate the data to 1H resolution
  #
  _year = dict["year"]
  _days = isleapyear(_year) ? 366 : 365
  @inline nt_to_1h(label::String) = (
    _dat_in = dict[label];
    @assert length(_dat_in) in [366, 365, 53, 52, 46, 12, 1] "Dataset length not supported";

    if length(_dat_in) == 1
      _dat_1d = repeat([_dat_in;]; inner=_days)
    elseif length(_dat_in) == 12
      _dat_1d = [([repeat(_dat_in[_m:_m], month_days(_year, _m)) for _m in 1:12]...)...]
    elseif length(_dat_in) == 46
      _dat_1d = repeat(_dat_in; inner=8)[1:_days]
    elseif length(_dat_in) in [52, 53]
      _dat_1d = repeat([_dat_in; _dat_in[end]]; inner=7)[1:_days]
    elseif length(_dat_in) in [365, 366]
      _dat_1d = [_dat_in; _dat_in[end]][1:_days]
    end;

    return repeat(_dat_1d; inner=24)
  )
  _df_in[!, "CO2"] .= nt_to_1h("co2_concentration")
  _df_in[!, "Chlorophyll"] .= nt_to_1h("chlorophyll")
  _df_in[!, "CI"] .= nt_to_1h("clumping_index")
  _df_in[!, "LAI"] .= nt_to_1h("leaf_area_index")
  _df_in[!, "Vcmax"] .= nt_to_1h("vcmax")

  # add the fields to store outputs
  for _label in DF_VARIABLES
    _df_in[!, _label] .= 0.0
  end

  return _df_in
end
