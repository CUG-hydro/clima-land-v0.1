using Revise
includet("main_pkgs.jl")
using RTableTools
# using NCDatasets

## 准备输入数据
f = "$(@__DIR__)/debug.nc"
fout = "$(@__DIR__)/debug.output.nc"

df_in = read_nc(f)
names(df_in)

@time dict = load("$(@__DIR__)/debug.jld2");
@time df = prepare_wd(dict, f); # vars: 27 -> 43


@time spac = prepare_spac(dict);
@time df_out = run_model!(spac, wddf; use_sif=false);
fwrite(df_out, "out-v0.csv")


FT=Float64

## 运行模型
_in_rad_bak = deepcopy(spac.in_rad)
_in_dir = _in_rad_bak.E_direct' * spac.wl_set.dWL / 1000
_in_dif = _in_rad_bak.E_diffuse' * spac.wl_set.dWL / 1000
_deepcopies = [_in_rad_bak, _in_dir, _in_dif]
_beta_g = BetaGLinearPsoil{FT}()

# set up memory
_spac_mem = SPACMemory{FT}()

_dfr = df[1, :]
prescribe_parameters!(spac, _dfr, _spac_mem, _deepcopies)
run_time_step!(spac, _dfr, _beta_g)

# # iterate through the time steps
# @showprogress for _dfr in eachrow(df)
#   prescribe_parameters!(spac, _dfr, _spac_mem, _deepcopies)
#   run_time_step!(spac, _dfr, _beta_g)
# end


# # save simulation results to hard drive
# save_nc!(nc_out, df[:, DF_VARIABLES])
df[:, DF_VARIABLES]
# return nothing
# @time df_out = run_model!(spac, wddf);
# @profview df_out = run_model!(spac, wddf);
