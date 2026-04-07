---
name: tadf-screening
description: TADF/蓝光分子高通量筛选技能。用于“先生成结构再筛选”的完整流程编排：10k 结构生成、xTB 预筛、TDDFT-xTB 波长窗口过滤、Gaussian TDDFT 精筛。触发场景：用户要求批量筛选、蓝光窗口过滤、远程 Marcus 筛选、阶段化漏斗筛选、进度巡检与异常恢复。
---

# TADF Screening (Standardized)

## 1) 强制前置：Workflow Initializer

在任何生成/筛选循环开始前，先运行：

```bash
python tadf-screening/scripts/screening_workflow_initializer.py
```

作用：
- 硬件/集群分层（`local_basic` / `local_gpu` / `remote_cluster`）
- 发光目标定义（波长范围、发光类型、谱宽）
- 数据库边界（DeepChem/PubChem + 自定义 `.csv/.smi`）
- TDDFT 引擎发现与选择（gaussian/orca/qchem/pyscf）
- 金字塔筛选参数（RDKit拼接+MMFF94+xTB稳定性 → xTB/sTDDFT+Stokes → 全 TDDFT）

> 关键化学约束：发射谱不能只用 S0 垂直激发。必须先做激发态几何优化（S1 或 T1）再算发射。

---

## 2) 标准流程（必须按顺序）

1. **结构生成（目标 10k）**
   - 输入：片段库/筛选计划
   - 输出：`structures/blue_plan_10k_xyz/*.xyz`

2. **xTB 预筛（Marcus 优先）**
   - 输入：`manifest.csv(idx,name,xyz_path)`
   - 输出：`xtb_progress.csv`, `xtb_state.json`

3. **TDDFT-xTB 波长窗口过滤（可选但推荐）**
   - 输入：xTB `status=ok` 候选
   - 输出：`tddft_xtb_results.csv`, `tddft_xtb_blue_window.csv`

4. **Gaussian TDDFT 精筛（终筛）**
   - 输入：第 3 步通过者
   - 协议：`S0 opt -> S1/T1 opt -> vertical emission`

---

## 3) 阶段门控（强制）

- 若 xTB `ok=0`：先排查环境/输入，禁止进入下一阶段。
- 若缺 `xtb4stda/stda`：记录 `tool_missing`，暂停 TDDFT-xTB 阶段。
- 若远程不可达：先恢复 SSH/路径，再重启筛选。

---

## 4) 执行位置标记（强制）

每次筛选开始前，必须在 `BATCH_STATUS.md` 记录：

- `xtb_execution_target: local|remote`
- `remote_host`（remote 时必填）
- `remote_workdir`（remote 时必填）
- `started_at`

未标记不得启动筛选。

---

## 5) 进度巡检（强制）

- 触发本 skill 后立即巡检一次。
- 任务进行中每 10 分钟至少巡检一次：
  - 进程是否存活
  - `xyz` 数量是否增长
  - `ok/fail/skip` 是否推进
- 异常时立刻恢复（重启/切换稳态脚本）并汇报。

---

## 6) 标准脚本入口

- 结构生成：
  - `tadf-screening/scripts/generate_10k_structures_fast.py`
- xTB 批筛：
  - `tadf-screening/scripts/run_xtb_batch_manifest.py`
- TDDFT-xTB 过滤：
  - `tadf-screening/scripts/run_tddft_xtb_filter.py`
- 前置初始化器：
  - `tadf-screening/scripts/screening_workflow_initializer.py`
- 拼接工具：
  - `tadf-screening/tools/smiles_assembler.py`
  - `tadf-screening/tools/molzip_assembler.py`

---

## 7) I/O 契约（最小）

### xTB 输入
```csv
idx,name,xyz_path
1,PXZ-TRZ,/abs/path/PXZ-TRZ.xyz
```

### xTB 输出
```csv
idx,name,xyz_path,status,detail,total_energy_eh,homo_lumo_gap_ev,normal_termination
```

### TDDFT-xTB 输出
```csv
idx,name,status,s1_ev,lambda_nm,delta_nm,pass_blue_window
```

---

## 8) 默认蓝光参数

- 目标波长窗口：`470 ± 30 nm`
- 窄带优先：FWHM 约 `0.15 eV`
- 广谱：FWHM 约 `0.30 eV`

如用户未指定，按以上默认执行。