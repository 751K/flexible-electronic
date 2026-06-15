# OTFT ECG 模拟前端（AFE）

基于 AT_4000TG PMOS-OTFT 工艺的**全差分 ECG 模拟前端**，包含 Cadence Virtuoso 电路与版图设计，以及一套与 Spectre 校准的 Python 求解器（DC + AC + Noise）。

## 项目结构

```
├── topology.py          # 电路拓扑定义（10 管全差分放大器）
├── pmos_tft_model.py    # AT_4000TG 器件模型（Verilog-A 等价 Python 实现）
├── ac_mna.py            # MNA 小信号 stamp 原语
├── ac_solver.py         # DC 工作点 + AC 小信号求解器
├── noise_solver.py      # 噪声求解器（输出噪声 + 输入参考噪声 IRN）
│
├── FGB/                 # 单元库：原理图 + 版图 + 符号 + 提取后版图
├── FGB_layout/          # 单元库（版图版本）
├── top/                 # 顶层原理图 + 符号
├── top_tb/              # 测试台（原理图 + Maestro 仿真配置）
└── Two/                 # 单元库（两管结构）
```

## 能做什么

- **Python 求解器**：在 Python 中直接运行 DC、AC、Noise 分析，结果与 Cadence Spectre 24.1 校准（典型点 <0.5%，最坏点 gain 差 0.01 dB、IRN 差 2%）。支持工艺角、逐管失配 Monte Carlo 扫描。
- **Cadence Virtuoso**：完整的电路原理图、版图（含提取后版图）和测试台，可直接在 Virtuoso 中打开仿真。

## 快速开始

**依赖**：Python 3 + NumPy + SciPy

```bash
pip install numpy scipy
```

**运行 DC + AC 分析**：

```bash
python ac_solver.py
```

**运行噪声分析**（含 IRN）：

```bash
python noise_solver.py
```

**在代码中使用**：

```python
import numpy as np
from noise_solver import noise_analysis, band_rms

sizes = {
    "M6": (2264, 78), "M7": (61365, 61), "M8": (61365, 61),
    "M9": (3175, 468), "M10": (3175, 468), "M11": (465, 66),
    "M12": (894, 85), "M13": (894, 85), "M14": (5224, 46), "M15": (5224, 46),
}
bias = {"VDD": 40.0, "VCM": 30.65, "VB": 9.84, "VC": 16.0}
freqs = np.logspace(-2, 4, 121)

r = noise_analysis(sizes, bias, freqs)
irn = band_rms(freqs, r["irn_psd"], 0.05, 100) * 1e6
print(f"IRN = {irn:.1f} µVrms")

# 指定叉指数（影响栅电容 → BW 和噪声）
nf = {"M7": 120, "M8": 120, "M9": 8, "M10": 8}
r = noise_analysis(sizes, bias, freqs, nf=nf)
```

## 器件模型与求解器文件关系

```
topology.py        ← 电路拓扑（无内部依赖）
pmos_tft_model.py  ← 器件模型（无内部依赖）
ac_mna.py          ← MNA stamp 原语（无内部依赖）
ac_solver.py       ← DC + AC（依赖以上三个）
noise_solver.py    ← Noise（依赖以上全部）
```

## 设计指标

| 指标 | 值 |
|------|-----|
| 增益 | 22.9 dB |
| 带宽 | 549 Hz |
| IRN（0.05–100 Hz） | 37 µVrms |
| 面积 | 11.3 mm² |
| 工艺角 × MC 良率 | 3 角 × 500 次，100% |

## 校准状态

与 **Cadence Spectre 24.1.0.078** 对照：
- 典型 / 工艺角：gain ±0.01 dB，BW 同点，IRN ≤ 数%
- 逐管失配 MC：均值与 σ 均与 Cadence 一致
- 8 个随机点全范围：最坏 gain 0.01 dB / BW 0% / IRN 2%

## 文档

技术细节（器件模型、求解器算法、校准方法等）见[技术文档](docs/README.md)。

## License

见 [LICENSE](LICENSE)。
