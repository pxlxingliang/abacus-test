# abacustest collectdata 提取指标参考

**版本**: abacustest v0.4.57  
**最后更新**: 2026-03-20

本文档列出 `abacustest collectdata` 命令支持的所有提取指标。使用 `-t` 参数指定数据类型：

| 值 | 数据类型 | 说明 |
|----|----------|------|
| `0` 或省略 | `abacus` | ABACUS 计算结果 |
| `1` | `qe` | Quantum ESPRESSO 计算结果 |
| `2` | `vasp` | VASP 计算结果 |

---

## ABACUS 指标 (`-t 0` 或 `--type abacus`)

### 基本信息

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `version` | ABACUS 版本 | 字符串 |
| `ncore` | MPI 进程数 | 整数 |
| `omp_num` | OMP 线程数 | 整数 |
| `normal_end` | 是否正常结束 | 布尔值 |
| `INPUT` | INPUT 文件参数 | 字典 |
| `kpt` | K 点设置 | 列表 |
| `fft_grid` | 电荷/势的 FFT 网格 | 列表 |
| `nbase` | LCAO 基组数量 | 整数 |
| `noccu_band` | 占据能带数 | 整数 |

### 系统参数

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `nbands` | 能带总数 | 整数 |
| `nkstot` | K 点总数 | 整数 |
| `ibzk` | 不可约 K 点数 | 整数 |
| `natom` | 原子总数 | 整数 |
| `nelec` | 电子总数 | 浮点数 |
| `nelec_dict` | 各元素电子数 | 字典 |
| `point_group` | 点群 | 字符串 |
| `point_group_in_space_group` | 空间群中的点群 | 字符串 |

### 磁性

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `total_mag` | 总磁矩 | Bohr mag/cell |
| `absolute_mag` | 绝对磁矩 | Bohr mag/cell |
| `total_mags` | 每个离子步的总磁矩 | 列表 |
| `absolute_mags` | 每个离子步的绝对磁矩 | 列表 |
| `atom_mag` | 每个原子的磁矩（最后离子步） | 列表 |
| `atom_mags` | 每个离子步的原子磁矩 | 列表 |
| `atom_orb_mag` | 原子轨道磁矩（最后离子步） | 列表 |
| `atom_orb_mags` | 每个离子步的原子轨道磁矩 | 列表 |
| `atom_elec` | 原子电子数（最后离子步） | 列表 |
| `atom_elecs` | 每个离子步的原子电子数 | 列表 |
| `atom_orb_elec` | 原子轨道电子数（最后离子步） | 列表 |
| `atom_orb_elecs` | 每个离子步的原子轨道电子数 | 列表 |
| `atom_mag_mul` | Mulliken 磁矩（最后离子步） | 列表 |
| `atom_mags_mul` | 每个离子步的 Mulliken 磁矩 | 列表 |
| `atom_elec_mul` | Mulliken 原子电子数 | 列表 |
| `atom_orb_elec_mul` | Mulliken 原子轨道电子数 | 列表 |
| `atom_elecs_mul` | 每个离子步的 Mulliken 原子电子数 | 列表 |
| `atom_orb_elecs_mul` | 每个离子步的 Mulliken 原子轨道电子数 | 列表 |

### 能量与收敛

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `converge` | SCF 是否收敛 | 布尔值 |
| `energy` | 总能量 | eV |
| `energy_ks` | Kohn-Sham 能量 | eV |
| `energies` | 每个离子步的总能量 | 列表 |
| `energy_per_atom` | 单原子能量 | eV/atom |
| `drho` | 每个 SCF 步的 drho | 列表 |
| `drho_last` | 最后 SCF 步的 drho | 浮点数 |
| `denergy` | 每个 SCF 步的 denergy | 列表 |
| `denergy_last` | 最后 SCF 步的 denergy | 浮点数 |
| `denergy_womix` | 每个 SCF 步的 denergy（未混合） | 列表 |
| `denergy_womix_last` | 最后 SCF 步的 denergy（未混合） | 浮点数 |

### 力与应力

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `force` | 原子力（最后离子步） | 列表 [3×natom] |
| `forces` | 每个离子步的原子力 | 列表 [nstep, 3×natom] |
| `stress` | 应力（最后离子步） | 列表 [9] |
| `stresses` | 每个离子步的应力 | 列表 [nstep, 9] |
| `virial` | 维里（最后离子步） | 列表 [9] |
| `virials` | 每个离子步的维里 | 列表 [nstep, 9] |
| `pressure` | 压强 | kbar |
| `pressures` | 每个离子步的压强 | 列表 |
| `largest_gradient` | 每个离子步的最大力梯度 | 列表 [eV/Å] |
| `largest_gradient_stress` | 每个离子步的最大应力梯度 | 列表 [kbar] |

### 结构与晶胞

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `volume` | 晶胞体积 | Å³ |
| `lattice_constant` | 晶格常数 [a,b,c,α,β,γ]（最后离子步） | 列表 |
| `lattice_constants` | 每个离子步的晶格常数 | 列表 |
| `cell` | 晶胞矢量（最后离子步） | 列表 [[3],[3],[3]] |
| `cells` | 每个离子步的晶胞矢量 | 列表 |
| `cell_init` | 初始晶胞矢量 | 列表 [[3],[3],[3]] |
| `coordinate` | 原子坐标（最后离子步） | 列表 [natom, 3] |
| `coordinates` | 每个离子步的原子坐标 | 列表 |
| `coordinate_init` | 初始原子坐标 | 列表 [natom, 3] |
| `element` | 元素名称列表 | 列表 |
| `label` | 原子标签列表 | 列表 |
| `element_list` | 同 element | 列表 |
| `atomlabel_list` | 同 label | 列表 |

### 能带与态密度

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `efermi` | 费米能级 | eV |
| `band` | 能带 | 列表 [nspin, nk, nband] |
| `band_weight` | 能带权重 | 列表 [nspin, nk, nband] |
| `band_plot` | 能带图文件名 | 字符串 |
| `e_bandgap` | 能隙（从 running_scf.log 输出） | eV |
| `band_gap` | 能隙 | eV |
| `k_coord` | K 点直接坐标 | 列表 |
| `dos` | 态密度（字典，含 energy 和 total） | 字典 |
| `pdos` | 投影态密度 | 字典 |

### 时间统计

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `total_time` | 总计算时间 | 秒 |
| `stress_time` | 应力计算时间 | 秒 |
| `force_time` | 力计算时间 | 秒 |
| `scf_time` | SCF 总时间 | 秒 |
| `scf_time_each_step` | 每个 SCF 步的时间 | 列表 |
| `step1_time` | 第一个 SCF 步的时间 | 秒 |
| `scf_steps` | SCF 迭代步数 | 整数 |

### 弛豫

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `relax_converge` | 弛豫是否收敛 | 布尔值 |
| `relax_steps` | 离子步总数 | 整数 |

### DeltaSpin

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `ds_lambda_step` | DeltaSpin 收敛步数 | 列表 |
| `ds_lambda_rms` | DeltaSpin RMS | 列表 |
| `ds_time` | DeltaSpin 内循环总时间 | 列表 |
| `ds_mag` | 原子磁矩 | 列表 |
| `ds_mags` | 每个离子步的原子磁矩 | 列表 |
| `ds_mag_force` | 原子磁力 | 列表 [eV/uB] |
| `ds_mag_forces` | 每个离子步的原子磁力 | 列表 |

### 内存

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `mem_vkb` | VNL::vkb 内存 | MB |
| `mem_psipw` | PsiPW 内存 | MB |

---

## Quantum ESPRESSO 指标 (`-t 1` 或 `--type qe`)

### 基本信息

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `normal_end` | 是否正常结束 | 布尔值 |
| `version` | QE 版本 | 字符串 |
| `ncore` | MPI 进程数 | 整数 |
| `natom` | 原子总数 | 整数 |

### 输入参数

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `ecutwfc` | 平面波截断能 | Ry |
| `ks_solver` | 对角化方法 | 字符串 |
| `mixing_type` | 电荷混合方法 | 字符串 |
| `mixing_beta` | 电荷混合系数 | 浮点数 |
| `scf_thr` | SCF 收敛阈值 | Ry |
| `smearing_method` | 展宽方法 | 字符串 |
| `smearing_sigma` | 展宽宽度 | Ry |

### 能量与收敛

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `converge` | SCF 是否收敛 | 布尔值 |
| `scf_steps` | SCF 迭代步数 | 整数 |
| `energy` | 总能量 | eV |
| `energies` | 每个离子步的总能量 | 列表 |
| `energy_per_atom` | 单原子能量 | eV/atom |

### 系统参数

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `nbands` | 能带数 | 整数 |
| `ibzk` | 不可约 K 点数 | 整数 |
| `nkstot` | K 点总数 | 整数 |
| `kpt` | K 点设置 | 列表 |
| `nelec` | 电子总数 | 浮点数 |

### 力与应力

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `force` | 原子力 | 列表 [eV/Å] |
| `forces` | 每个离子步的原子力 | 列表 |
| `stress` | 应力 | 列表 [kbar] |
| `stresses` | 每个离子步的应力 | 列表 |
| `virial` | 维里 | 列表 [eV] |
| `virials` | 每个离子步的维里 | 列表 |
| `pressure` | 压强 | kbar |
| `pressures` | 每个离子步的压强 | 列表 |

### 结构与磁性

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `cell` | 晶胞矢量 | 列表 [Å] |
| `volume` | 晶胞体积 | Å³ |
| `coord` | 原子坐标 | 列表 [Å] |
| `label` | 原子标签 | 列表 |
| `atomlabel_list` | 同 label | 列表 |
| `element` | 元素名称 | 列表 |
| `element_list` | 同 element | 列表 |
| `total_mag` | 总磁矩 | 浮点数 |
| `absolute_mag` | 绝对磁矩 | 浮点数 |
| `atom_mag` | 原子磁矩 | 列表 |

### 能带

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `band` | 能带 | 列表 |
| `efermi` | 费米能级 | eV |
| `band_gap` | 能隙 | eV |

### 收敛与时间

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `drho` | 每个 SCF 步的 drho | 列表 [Ry] |
| `drho_last` | 最后 SCF 步的 drho | 浮点数 |
| `denergy` | 每个 SCF 步的 denergy | 列表 [eV] |
| `denergy_last` | 最后 SCF 步的 denergy | 浮点数 |
| `total_time` | 总时间（WALL） | 秒 |
| `force_time` | 力计算时间 | 秒 |
| `stress_time` | 应力计算时间 | 秒 |
| `scf_time` | SCF 总时间 | 秒 |
| `scf_time_each_step` | 每个 SCF 步的时间 | 列表 |
| `scf_time_step1` | 第一个 SCF 步的时间 | 秒 |

### 弛豫

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `relax_converge` | 弛豫是否收敛 | 布尔值 |
| `relax_steps` | 离子步总数 | 整数 |

---

## VASP 指标 (`-t 2` 或 `--type vasp`)

### 基本信息

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `version` | VASP 版本 | 字符串 |
| `ncore` | MPI 进程数 | 整数 |
| `normal_end` | 是否正常结束 | 布尔值 |

### 输入参数

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `encut` | 截断能 | eV |
| `ismear` | 展宽方法 | 整数 |
| `sigma` | SIGMA 值 | eV |
| `nelm` | 最大 SCF 步数 | 整数 |
| `natom` | 原子总数 | 整数 |
| `nbands` | 能带数 | 整数 |
| `nelec` | 电子总数 | 浮点数 |
| `spin` | 自旋数 | 整数 |
| `volume` | 体积（最后离子步） | Å³ |

### DFT+U 参数

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `ldautype` | LDAUTYPE 值 | 整数 |
| `ldaul` | LDAUL 值（各元素） | 列表 |
| `ldauu` | LDAUU 值（各元素） | 列表 |
| `ldauj` | LDAUJ 值（各元素） | 列表 |

### K 点

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `kpt` | K 点设置 | 列表 |
| `nkstot` | K 点总数 | 整数 |
| `ibzk` | 不可约 K 点数 | 整数 |

### 能量与收敛

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `converge` | SCF 是否收敛 | 布尔值 |
| `scf_steps` | SCF 迭代步数 | 整数 |
| `energy` | 总能量（最后离子步） | eV |
| `energies` | 每个离子步的总能量 | 列表 |
| `energy_per_atom` | 单原子能量 | eV/atom |
| `denergy` | SCF 能量差 | eV |
| `denergy_last` | 最后 SCF 步的能量差 | eV |

### 力与应力

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `force` | 原子力 | 列表 [eV/Å] |
| `forces` | 每个离子步的原子力 | 列表 |
| `stress` | 应力 | 列表 [kbar] |
| `stresses` | 每个离子步的应力 | 列表 |
| `virial` | 维里 | 列表 [eV] |
| `virials` | 每个离子步的维里 | 列表 |
| `pressure` | 压强 | kbar |
| `pressures` | 每个离子步的压强 | 列表 |

### 结构与元素

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `cell` | 晶胞矢量（最后离子步） | 列表 [[3],[3],[3]] |
| `cells` | 每个离子步的晶胞矢量 | 列表 |
| `atom_name` | 元素名称 | 列表 |
| `atom_type` | 元素类型名称 | 列表 |
| `element` | 元素名称列表 | 列表 |
| `element_list` | 同 element | 列表 |
| `label` | 原子标签列表 | 列表 |
| `atomlabel_list` | 同 label | 列表 |

### 磁性与能带

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `total_mag` | 总磁矩 | 浮点数 |
| `absolute_mag` | 绝对磁矩 | 浮点数 |
| `atom_mag` | 原子磁矩 | 列表 |
| `efermi` | 费米能级 | eV |
| `band` | 能带 | 列表 [spin×kpoint×band] |
| `band_gap` | 能隙 | eV |
| `band_plot` | 能带图文件名 | 字符串 |

### 对称性

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `point_group` | 点群 | 字符串 |
| `point_group_in_space_group` | 空间群中的点群 | 字符串 |

### 时间与弛豫

| 指标名 | 说明 | 单位/类型 |
|--------|------|-----------|
| `total_time` | 总 CPU 时间 | 秒 |
| `scf_time` | SCF 总时间 | 秒 |
| `stress_time` | 应力计算时间 | 秒 |
| `relax_steps` | 离子步数 | 整数 |
| `relax_converge` | 弛豫是否收敛 | 布尔值 |

---

## 使用示例

### 提取 ABACUS 结果

```bash
# 提取所有默认指标
abacustest collectdata -j job_dir -o results.json

# 提取特定指标
abacustest collectdata -j job_dir -t 0 -p energy,force,stress -o results.json

# 提取多个任务
abacustest collectdata -j job1 job2 job3 -o all_results.json
```

### 提取 QE 结果

```bash
abacustest collectdata -j qe_job -t 1 -o qe_results.json
```

### 提取 VASP 结果

```bash
abacustest collectdata -j vasp_job -t 2 -o vasp_results.json
```

---

## 注意事项

1. **delta_ 前缀指标**: 以 `delta_` 开头的指标需要参考值 JSON 文件，通过 `--ref` 指定
2. **模块限制**: 默认只使用主模块（abacus/qe/vasp），可通过 `--modules` 使用其他模块
3. **单位**: 注意不同 DFT 软件的默认单位可能不同
4. **列表维度**: 注意 `force`、`stress` 等列表的维度说明

---

*本文档由 abacustest collectdata --outparam 命令输出整理*
