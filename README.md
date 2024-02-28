# CliMA Land v0.1

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/clima-land-v0.1/dev)
[![CI](https://github.com/CUG-hydro/clima-land-v0.1/actions/workflows/JuliaStable.yml/badge.svg)](https://github.com/CUG-hydro/clima-land-v0.1/actions/workflows/JuliaStable.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/clima-land-v0.1/branch/test/graph/badge.svg)](https://codecov.io/gh/CUG-hydro/clima-land-v0.1)


The land component of CliMA ESM is moving towards broadband RT and big leaf modeling, which does not align with my research interest to make the land model more complicated but powerful, such as using remote sensing data to calibrate the model. Therefore, I create this new repository to continue my development.

This respository clones the final code of CliMA Land v0.1 rather than the commit history. Please refer to [CliMA Land](https://github.com/CliMA/Land) branch v0.1 for the original land model and other branches for its future developments.


## About
At the moment this repository is created, there are two land models within CliMA: [CliMA Land](https://github.com/CliMA/Land) and [ClimaLSM.jl](https://github.com/CliMA/ClimaLSM.jl). The two models have different aims that the former prefers complex models and use more field- and satellite-based observations, while the latter prefers simple models and runs faster to serve the Atmosphere and Ocean models developed within CliMA. There might be future plans to merge the two repository, but it is not clear when and how it will happen. Therefore, to not complicate the problem, I am moving my development to this new repository. Hopefully in the future, CliMA software team will have some time to help the migration of the complex model we have developed into CliMA.

Another note is that this v0.1 CliMA Land is relatively outdated, given that my time has been spent on v0.2 of CliMA Land, which has a very different framework (the main branch of CliMA Land when creating this repository). And for the same reason of creating this repository, future developments are made on [Emerald](https://github.com/Yujie-W/Emerald). Therefore, this repository is mainly for bug fixes and sometimes adding new features before Emerald is ready for global simulations.


## Examples

### Run CliMA Land for a single site (v0.1)

1. Download all the files including the source code in folder examples into a folder,
```shell
$ git clone https://github.com/silicormosia/clima-land-v0.1
$ cd clima-land-v0.1
$ cd examples
```

2. initialize the Julia environment
```shell
$ julia --project -e "using Pkg; Pkg.instantiate();"
```

3. Make sure you have these files in the folder before heading to next step
   - `Project.toml`
   - `Manifest.toml`
   - `example.jl`
   - `debug.jld2`
   - `debug.nc`

4. run the model for the site
```shell
$ julia --project example.jl
```

5. you should get a file named `debug.output.nc` after a few minutes

You will need to edit the functions we provided if you need to change more parameters, or output more results. Feel free to contact us through email or Github Issues (if this tutorial does not work).


## References

Please cite the following when you use the CliMA Land (v0.1)

### General model description
Y. Wang, P. Köhler, L. He, R. K. Braghiere, R. Doughty, J. Wood, C. Frankenberg. 2021.
Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1).
Geoscientific Model Development. 14(11): 6741-6763.
[DOI](https://doi.org/10.5194/gmd-14-6741-2021)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing.pdf)
[SI](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2021testing-si.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
@article{wang2021testing,
    author = {Wang, Y. and K{\"o}hler, P. and He, L. and Doughty, R. and Braghiere, R. K. and Wood, J. D. and Frankenberg, C.},
    year = {2021},
    title = {Testing stomatal models at the stand level in deciduous angiosperm and evergreen gymnosperm forests using CliMA Land (v0.1)},
    journal = {Geoscientific Model Development},
    volume = {14},
    number = {11},
    pages = {6741--6763}
}
```

### Clumping index implementation
R. K. Braghiere, Y. Wang, R. Doughty, D. Souza, T. Magney, J. Widlowski, M. Longo, A. Bloom, J. Worden, P. Gentine, and C. Frankenberg. 2021.
Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model.
Remote Sensing of Environment. 261: 112497.
[DOI](https://doi.org/10.1016/j.rse.2021.112497)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting.pdf)
[SI](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/braghiere2021accounting-si.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
@article{braghiere2021accounting,
    author = {Braghiere, Renato K and Wang, Yujie and Doughty, Russell and Sousa, Daniel and Magney, Troy and Widlowski, Jean-Luc and Longo, Marcos and Bloom, A Anthony and Worden, John and Gentine, Pierre and Frankenberg, Christian},
    year = {2021},
    title = {Accounting for canopy structure improves hyperspectral radiative transfer and sun-induced chlorophyll fluorescence representations in a new generation Earth System model},
    journal = {Remote Sensing of Environment},
    volume = {261},
    pages = {112497}
}
```

### Canopy complexity
Y. Wang, C. Frankenberg. 2022.
On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes.
Biogeosciences. 19(1): 29-45.
[DOI](https://doi.org/10.5194/bg-19-29-2022)
[PDF](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG/raw/master/publications/wang2022impact.pdf)
[CODE](https://github.com/Yujie-WANG/Published-Codes-Yujie-WANG)

```
 @article{wang2022impact,
 	author = {Wang, Y. and Frankenberg, C.},
 	year = {2022},
 	title = {On the impact of canopy model complexity on simulated carbon, water, and solar-induced chlorophyll fluorescence fluxes},
 	journal = {Biogeosciences},
 	volume = {19},
 	number = {1},
 	pages = {29--45}
 }
 ```
