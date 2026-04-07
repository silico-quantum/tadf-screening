---
name: tadf-screening
version: 1.0.0
description: TADF 发射体高通量筛选工具。支持从 SMILES 批量生成 Gaussian 输入，自动化计算和结果提取。
homepage: https://github.com/your-repo/tadf-screening
metadata:
  category: computational_chemistry
  tags: [TADF, screening, TD-DFT, STO-3G, automation]
---

# TADF 高通量筛选技能

## 概述

用于 TADF (Thermally Activated Delayed Fluorescence) 发射体的高通量虚拟筛选。支持从 SMILES 字符串到 Gaussian 计算的完整流程。

## 核心功能

1. **批量输入文件生成**: SMILES → Gaussian input (.gjf)
2. **多种计算类型**: S0 优化、TD-DFT、T1 优化
3. **自动化作业提交**: Slurm 脚本生成
4. **结果提取**: ΔEST、λPL、f12 等关键性质

## 计算流程

### 完整协议 (Recommended)

```
1. SMILES → 3D 结构 (RDKit)
2. S0 几何优化 (B3LYP/STO-3G)
3. TD-DFT 垂直激发 (B3LYP/STO-3G)
4. T1 几何优化 (TD-DFT/STO-3G, root=1)
5. 估算 Stokes shift
6. 计算 ΔEST = E(S1) - E(T1)
```

### 快速协议 (High-throughput)

```
1. SMILES → 3D 结构
2. S0 几何优化
3. TD-DFT 垂直激发
4. 估算 ΔEST (vertical)
```

## 强制巡检规则（立即生效）

> 只要触发本 skill 并处于筛选/生成执行阶段，必须执行该规则。

1. **立即巡检**：触发 skill 后立刻做一次进程与产出检查（不等待）。
2. **10 分钟周期巡检**：在任务未完成前，每 10 分钟至少检查一次：
   - 运行会话是否存在（`process list`）
   - 结构数量是否增长（如 `blue_plan_10k_xyz/*.xyz`）
   - manifest 的 `ok/fail/skip` 是否推进
3. **异常处置**：若发现进程停止、卡死或无增长，必须立刻执行恢复动作（重启任务/切换稳态脚本），并在回复中说明。
4. **状态汇报**：每次巡检后都给出最小状态集：`当前数量 / 目标数量`、`进程状态`、`是否采取恢复动作`。
5. **完成条件**：达到目标数量（如 10k）或用户明确停止后，结束周期巡检。

## 标准筛选输入/输出（xTB 预筛）

### 标准输入

1. `manifest.csv`（批次清单，最少字段）
```csv
idx,name,xyz_path
1,PXZ-COTRZ__direct__A-linker-D__secondary,/abs/path/to/mol.xyz
```

2. `xyz_path` 指向有效 3D `.xyz` 文件（由上游结构生成提供）

### 标准执行

```bash
python tadf-screening/scripts/run_xtb_batch_manifest.py \
  --manifest <batch_manifest.csv> \
  --batch-dir <batch_dir> \
  --xtb-bin /opt/homebrew/bin/xtb \
  --max-items 5000 \
  --timeout 120 \
  --gfn 2
```

> 当前环境标准模式为 **single-point**：`xtb --gfn 2 --sp`。
> 原因：本机 xTB 在 `--opt` 下存在 Fortran runtime 崩溃，导致系统性失败。

### 执行位置标记（强制）

在**每次筛选开始前**，必须先声明 xTB 执行位置，并写入批次状态：

- `xtb_execution_target: local`（本地执行）
- `xtb_execution_target: remote`（远程服务器执行，如 marcus/marcus2）

同时在 `skills/tadf-screening/BATCH_STATUS.md` 的当前批次条目中至少记录：

- `xtb_execution_target`
- `remote_host`（当 target=remote 时必填，例如 `marcus2`）
- `remote_workdir`（当 target=remote 时必填）
- `started_at`

未完成该标记，不得启动筛选进程。

### 标准输出

1. `xtb_progress.csv`（逐分子结果）
```csv
idx,name,xyz_path,status,detail,total_energy_eh,homo_lumo_gap_ev,normal_termination
1,PXZ-TRZ,...,ok,ok,-59.170853,0.747637,True
```

2. `xtb_state.json`（断点续跑状态）
```json
{"next_index": 120, "done": 120, "ok": 118, "fail": 2}
```

3. `xtb_results/<name>/xtb.out`（原始日志）

### 通过判据（标准）

- `returncode == 0`
- `xtb.out` 中包含 `normal termination of xtb`

## 蓝光筛选标准流程（结构 → xTB → TDDFT-xTB → Gaussian）

按以下顺序执行，不跳步：

1. **生成结构池（目标 10k）**
   - 输出：`blue_plan_10k_xyz/*.xyz`
   - 规则：先 SMILES/拼接合法性，再 3D 结构生成

2. **xTB 预筛（几何与电子粗筛）**
   - 输入：`manifest.csv (idx,name,xyz_path)`
   - 输出：`xtb_progress.csv`, `xtb_state.json`
   - 只保留 `status=ok` 候选进入下一步

3. **TDDFT-xTB 波长窗口过滤（蓝光）**
   - 输入：xTB 的 `status=ok` 候选
   - 目标：按目标发光波长窗口过滤（默认 `470 ± 30 nm`）
   - 输出：`tddft_xtb_results.csv`, `tddft_xtb_blue_window.csv`

4. **Gaussian TDDFT 精筛**
   - 仅对 `tddft_xtb_blue_window.csv` 通过者提交 Gaussian
   - 评估 `ΔEST`、`f12`、发光波长等最终指标

### 阶段门控（强制）

- 若第 2 步 xTB `ok=0`：必须先排查输入/环境，禁止直接进入第 3/4 步。
- 若第 3 步工具缺失（`xtb4stda/stda`）：必须记录 `tool_missing` 并暂停该步，不得伪造结果。
- 每步开始前必须更新 `BATCH_STATUS.md`：批次、执行位置（local/remote）、输入、输出、启动时间。

## 新增阶段：TDDFT-xTB 波长窗口过滤（蓝光）

在 xTB 预筛后增加一层激发态过滤，尽早排除与目标发光波长偏差太大的分子。

### 标准输入

- `xtb_progress.csv`（来自 xTB 预筛，至少包含 `status,name,xyz_path`）
- 目标波长参数：`target_nm`（默认 470）与 `window_nm`（默认 ±30）

### 标准执行

```bash
python tadf-screening/scripts/run_tddft_xtb_filter.py \
  --xtb-progress <xtb_progress.csv> \
  --outdir <tddft_xtb_dir> \
  --target-nm 470 \
  --window-nm 30 \
  --max-items 500
```

> 说明：TDDFT-xTB 通常依赖 `xtb4stda + stda`。若环境缺失工具，脚本会显式输出 `tool_missing`，不静默失败。

### 标准输出

- `tddft_xtb_results.csv`：全量结果（含失败原因、S1 eV、估算波长）
- `tddft_xtb_blue_window.csv`：通过蓝光窗口过滤的候选

### 过滤规则

- 从输出中解析第一激发能 `S1(eV)`
- 估算波长：`lambda_nm = 1239.84 / S1(eV)`
- 保留：`abs(lambda_nm - target_nm) <= window_nm`

## 快速开始

### 1. 准备 SMILES 列表

创建 `molecules.csv`:
```csv
name,smiles
2CzPN,c1ccc2c(c1)ccc1cc3ccccc3cc-2n1
DMAC-TRZ,CN(C)c1ccc(Cn2cnc3cc(Sc4ccccc4)c(C)cc3c2=O)cc1
```

### 2. 生成输入文件

```bash
python3 scripts/generate_inputs.py molecules.csv --method b3lyp --basis sto-3g --type full
```

### 3. 提交作业

```bash
bash scripts/submit_jobs.sh tadf_jobs/
```

### 4. 提取结果

```bash
python3 scripts/extract_results.py tadf_jobs/ --output results.csv
```

## 计算类型

### opt_s0: S0 态几何优化

**输入模板**: `templates/opt_s0.gjf.template`

```gjf
%chk={name}_s0.chk
%mem=4GB
%nproc=4
# {method}/{basis} opt

{name} - S0 geometry optimization
0 1
{coordinates}
```

**用途**: 优化基态几何结构

### td_vertical: TD-DFT 垂直激发

**输入模板**: `templates/td_vertical.gjf.template`

```gjf
%chk={name}_td.chk
%mem=4GB
%nproc=4
# td=(nstates=10) {method}/{basis}

{name} - Vertical TD-DFT
0 1
{coordinates}
```

**用途**: 计算垂直激发能和振子强度

### opt_t1: T1 态几何优化

**输入模板**: `templates/opt_t1.gjf.template`

```gjf
%chk={name}_t1.chk
%mem=4GB
%nproc=4
# td=(nstates=10,root=1) {method}/{basis} opt

{name} - T1 geometry optimization
0 1
{coordinates}
```

**用途**: 优化三重态几何，用于计算绝热激发能

## 输出格式

### CSV 格式 (results.csv)

```csv
name,smiles,E_S0,E_S1,E_T1,delta_EST,eV,lambda_PL,nm,f12,oscillator_strength,calc_time,hours
2CzPN,c1ccc2c(c1)...,-1234.5678,-1234.2345,-1234.2456,0.0111,0.302,410.2,0.0234,2.5
```

## 关键性质定义

### ΔEST (Singlet-Triplet Gap)

- **定义**: ΔEST = E(S1) - E(T1)
- **单位**: eV (1 eV = 27.2114 Hartree)
- **TADF 标准**: ΔEST < 0.3 eV

### λPL (Photoluminescence Wavelength)

- **定义**: λPL = 1239.84 / E(S1 - S0)
- **单位**: nm
- **可见光范围**: 380-780 nm

### f12 (Oscillator Strength)

- **定义**: S0 → S1 跃迁的振子强度
- **意义**: 辐射跃迁速率的度量
- **典型 TADF**: f12 < 0.1 (CT 特征)

## 与 xTB 方法的对比

| 性质 | xTB (论文) | STO-3G (我们的) |
|------|-----------|----------------|
| ΔEST MAE | ~0.17 eV | 待验证 |
| λPL MAE | ~78.8 nm | 待验证 |
| 计算时间/分子 | ~1 min | ~5-30 min |
| 基组/方法 | Semi-empirical | Ab initio |

## 常见问题

### Q1: 计算时间太长?

**解决**:
- 减少激发态数量: `nstates=5` (默认 10)
- 使用更快的泛函: `pbe0` (vs `b3lyp`)
- 考虑 xTB 预筛选

### Q2: ΔEST 为负?

**可能原因**:
- 基态不是 S0 (自旋污染)
- 几何未收敛
- 需要检查波函数稳定性

**解决**:
- 添加 `stable=opt` 关键词
- 检查 <S²> 值

## 相关 Skills

- **gaussian**: Gaussian 计算的详细指南
- **multiwfn**: 波函数分析和轨道可视化
- **molecular-sampler**: 分子采样和构象生成

## 参考文献

1. 论文原文: arXiv:2511.00922v1
2. TADF 理论: Adachi et al., Nature 2012
3. TD-DFT 方法: Runge & Gross, PRL 1984

---

**维护者**: Silico (硅灵)
**创建时间**: 2026-04-02
