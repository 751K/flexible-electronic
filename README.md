# AFE 核心求解器栈（已对 Cadence 校准）

OTFT ECG AFE 的最小自洽 **DC + AC + Noise** 求解栈，对 Cadence Spectre 24.1 校准：
典型点 <0.5%，工艺角逐项匹配，逐管失配 MC 的 σ 与 Cadence 一致，随机点全范围
（含多稳态/退化/超轨边缘点）最坏偏差 **gain 0.01 dB / IRN 2%**。

**5 个文件**，按 import 关系闭合，直接 `python noise_solver.py` 即可跑（依赖 numpy + scipy）：

```
topology.py        ← 电路拓扑「单一真源」（无内部依赖）
pmos_tft_model.py  ← 器件模型（无内部依赖）
ac_mna.py          ← MNA stamp 原语（无内部依赖）
ac_solver.py       ← DC + AC   (import: topology, ac_mna, pmos_tft_model)
noise_solver.py    ← Noise     (import: ac_solver, topology, ac_mna, pmos_tft_model)
```

## 文件职责

### 1. `pmos_tft_model.py` — PMOS-OTFT 器件模型
AT_4000TG 的 Verilog-A 等价 Python 模型。`get_Idc` 解端子电流（内部含接触网络节点 s1/d1），
`get_noise_psd` 给漏电流噪声 PSD（Hooge flicker + thermal），`get_capacitances` 给 Cgs/Cgd，
`g_area` 给绘制面积。构造参数含**工艺角/失配旋钮** `pvt0/mvt0/pbeta0/mbeta0`。
> AC/噪声一律用**端子** gm/gds（对 `get_Idc` 做有限差分），不是沟道 gm。

### 2. `topology.py` — 电路拓扑「单一真源」★
`Topology` 类：一张器件表 `(name, drain, gate, source)`（按节点名）+ solved 节点表 + rails。
**DC KCL、逐管偏置映射、AC/噪声端子表全部从它派生**，不再在求解器里手写（消除了重复与抄写 bug）。
默认实例 `AFE_TOPO`（10 管全差分：M6 尾电流、M7/8 输入对、M9/10 输出级、M11–M15 交叉耦合正反馈电平移位）。
求解器经 `topo=AFE_TOPO` 形参接收，可换不同拓扑。

### 3. `ac_mna.py` — MNA stamp 原语
`_stamp_adm` / `_stamp_vccs` / `_stamp_mos`，整个栈的小信号引擎。

### 4. `ac_solver.py` — DC + AC
- `ac_solve(sizes, bias, freqs, corner=None, x0_guess=None, topo=AFE_TOPO)`
- DC：`scipy.fsolve` 解 6 节点（KCL 由 topology 派生）；返回增益曲线 / −3dB BW / **全节点 DC op** / 端子小信号 `ss`。
- **工艺角 + 逐管失配**：`corner` 可为
  - `None` → 标称；
  - 扁平 dict `{'pvt0':..,'pbeta0':..}` → 全局工艺角（所有管一致）；
  - 每管 map `{'M7':{...},...}` → 逐管失配（每管独立 mvt0/mbeta0）。
- **DC 鲁棒性三守护**（贴合 Spectre 的物理分支，避开多稳态假解）：
  ① **源步进连续法**主种子（模拟上电，追 Spectre 物理分支）；
  ② **对称守护**（无失配时若解对称破缺 latch，重解对称系统）；
  ③ **物理性守护**（节点超出 [0,VDD] 的非物理分支重解到轨内）。
  全部**只在冷启动/异常时触发**；带种子的 in-loop 扫描（MC/优化）走快路径，速度不受影响。
- `get_ss_params`：有限差分提取端子 gm/gds + 端子电容。

### 5. `noise_solver.py` — Noise
在与 AC **同一个** 6 节点 MNA 上做噪声传播（Spectre 同款方法）：每管漏电流噪声 PSD 注在其 drain/source，
读到差分输出 (vop−von) 的转移阻抗 Z_k，输出噪声 PSD = Σ|Z_k|²·S_id,k（不相关），再除以增益 |H| 折合到输入得 IRN。
偏置映射与端子表同样由 `topology` 派生（支持工艺角/失配透传）。

## 快速使用
```python
import numpy as np
from noise_solver import noise_analysis, band_rms
from ac_solver import ac_solve

sizes = {  # (W, L) in µm
    "M6": (2264,78), "M7": (61365,61), "M8": (61365,61),
    "M9": (3175,468), "M10": (3175,468), "M11": (465,66),
    "M12": (894,85), "M13": (894,85), "M14": (5224,46), "M15": (5224,46),
}
bias = {"VDD":40.0, "VCM":30.65, "VB":9.84, "VC":16.0}
freqs = np.logspace(-2, 4, 121)

ac  = ac_solve(sizes, bias, freqs)                  # 典型：gain/BW/DC op
r   = noise_analysis(sizes, bias, freqs)            # 噪声谱 + IRN
irn = band_rms(freqs, r["irn_psd"], 0.05, 100)*1e6  # µVrms

# 工艺角（slow）：
slow = {"pvt0": -3*0.0753, "pbeta0": -15*0.036}
r_slow = noise_analysis(sizes, bias, freqs, corner=slow)
```

## 校准状态（vs Cadence Spectre 24.1.0.078）
- 典型 / 工艺角：gain ±0.01 dB、BW 同点、IRN ≤ 数 %。
- 逐管失配 MC：均值与 σ 均与 Cadence 一致。
- 8 个随机点全范围对照：最坏 gain 0.01 dB / BW 0% / IRN 2%。
- 锁定的最终设计 11.3 mm²：22.9 dB / 549 Hz / 37 µV，三项达标，3 工艺角×500 MC 100% 良率。
