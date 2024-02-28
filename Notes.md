
# 现阶段，`SoilPlantAirContinuum`尚没有内容。

# 主要过程集中在`StomataModels`

```julia
gas_exchange!(spac.photo_set, _iPS, _iEN, GswDrive())
update_gsw!(spac, spac.stomata_model, _i_can, FT(120); β=_βm) # 求解gsw, 更新叶片导度
# prognostic_gsw!(spac.plant_ps[_i_can], spac.envirs[_i_can], spac.stomata_model, βm, FT(120))
gsw_control!(spac.photo_set, _iPS, _iEN)
```

> 三个需要学习

## `gas_exchange!`

计算光合速率


# TODO

- 土壤蒸发

- 冠层截留

- surface_fluxes_balance
  + 土壤温度

- 土壤水热过程

- 辐射模块：大叶模型、双叶模型
  + 叶片温度

- 地下水


# 学习部分

- 冠层辐射传输
