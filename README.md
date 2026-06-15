# current_core — 当前 AFE 设计的核心代码快照

这是 OTFT ECG AFE 当前**已对 Cadence 校准**的最小自洽代码栈（DC/AC/Noise 三项
误差 <0.5%，并在 6 组尺寸上交叉验证 gain ±0.01 dB / BW 同点 / IRN ≤1.6%）。

4 个文件按 import 关系闭合，直接 `python noise_solver.py` 即可跑（依赖 numpy + scipy）：

```
noise_solver.py ──> ac_solver.py ──> pmos_tft_model.py
              └──> ac_mna.py    ──> pmos_tft_model.py
```

## 三个组成部分 ↔ 文件对应

### 1. PMOS 模型
- **`pmos_tft_model.py`** — AT_4000TG PMOS-OTFT 的 Verilog-A 等价 Python 模型。
  内部含接触网络节点 (s1/d1)，`get_Idc` 解端子电流，`get_noise_psd` 给沟道漏电流
  噪声 PSD（Hooge flicker + thermal），`get_capacitances` 给 Cgs/Cgd。
  > 关键：AC/噪声用的是**端子** gm/gds（有限差分 `get_Idc`），不是沟道 gm。

### 2. 电路拓扑
拓扑没有独立文件，而是以**数据**内嵌在求解器里（10 管全差分：M6 尾电流、
M7/M8 输入对、M9/M10 输出级、M11–M15 交叉耦合正反馈电平移位）：
- **`ac_solver.py` 的 `residuals()`** — DC 工作点的 KCL 方程 = 直流拓扑（哪管接哪个节点）。
- **`ac_solver.py` 的 `devs` 列表** — 小信号 `(name, drain, gate, source)` 连接表 = 交流拓扑。
- **`noise_solver.py` 的 `devs` / `bpts`** — 同一拓扑 + 每管偏置映射（噪声注入用）。
- 节点编号统一：`VOP=0, VON=1, vfbp=2, vfbn=3, net20=4, net2=5`；负载 CL=5pF。

> 拓扑的文字版/原理图见 `../docs/step2_circuit_topology.md`；
> Cadence 网表见服务器 `~/afe_gt/tb.scs`（参数化：CurrentW/L→M6, InputW/L→M7/8,
> LoadW/L→M9/10, LevelW/L→M14/15, PW/PL→M12/13, cw/cl→M11）。

### 3. 求解器
- **`ac_mna.py`** — 干净的 MNA 小信号 stamp 原语：`_stamp_mos`（VCCS 跨导 + gds + Cgs/Cgd）、
  `_stamp_adm`（导纳）。整个栈的小信号引擎。
- **`ac_solver.py`** — 全电路 DC（scipy.fsolve 解 6 节点）+ AC（逐频点解 6 节点 MNA）。
  返回增益曲线、-3dB BW、DC 工作点、每管端子小信号参数 `ss`。
- **`noise_solver.py`** — 在 *同一个* AC MNA 上做噪声传播（Spectre 同款方法）：
  每管漏电流噪声 PSD 注在 drain/source，读到差分输出的转移阻抗 Z_k，
  输出 PSD = Σ|Z_k|²·S_id,k，再除以增益折算到输入得 IRN。
  > 噪声 MNA = AC MNA：噪声分析时输入无信号，M7/M8 栅是 AC 地，Y 矩阵与 AC 完全相同。

## 快速使用

```python
import numpy as np
from noise_solver import noise_analysis, band_rms
from ac_solver import ac_solve

sizes = {  # (W, L) in um
    "M6": (3000,150), "M7": (25000,150), "M8": (25000,150),
    "M9": (12000,500), "M10": (12000,500), "M11": (300,100),
    "M12": (500,80), "M13": (500,80), "M14": (2000,500), "M15": (2000,500),
}
bias  = {"VDD":40.0, "VCM":32.0, "VB":20.0, "VC":26.0}
freqs = np.logspace(-2, 4, 121)

ac = ac_solve(sizes, bias, freqs)              # gain / BW / DC op
r  = noise_analysis(sizes, bias, freqs)        # 噪声谱 + IRN
irn = band_rms(freqs, r["irn_psd"], 0.05, 100) # 等效输入噪声 [Vrms]
```

## 当前结论（截至 2026-05-29）
当前尺寸下：gain 19.96 dB ✅，但 **BW 50 Hz（规格 ≥100Hz，欠 2×）❌**、
**IRN 209 µV（规格 ≤44.5µV，欠 4.7×）❌**。且降噪要增大 WL → 节点电容 ∝ WL →
带宽进一步下降，二者在该拓扑下**直接冲突**，纯尺寸无可行解。详见 `../docs/design_state.md`。
