# flexible-electronic

## 目标指标

使用 printed Organic Thin-Film Transistors（OTFTs）设计一个用于 ECG sensing 的模拟前端电路（AFE）。ECG 信号由两片 flexible dry electrodes 进行差分采集，并通过 AC-coupled AFE 进行读出和放大。

### 应用需求

| 项目 | 目标 |
|---|---:|
| 最低信噪比 SNR | ≥ 18 dB |
| ECG 检测带宽 | 0.05 Hz – 100 Hz |
| SNR 计算方式 | 基于采样后的信号计算 |

### 电路性能指标

| 项目 | 目标 |
|---|---:|
| 差分输入信号幅度 | 0.5 mVpeak – 5 mVpeak |
| 差分 AFE 增益 | ≥ 10 V/V |
| 差分 AFE 增益，dB 表示 | ≥ 20 dB |
| AFE 3-dB 带宽 | 至少覆盖 0.05 Hz – 100 Hz |
| 输出采样频率 | ≥ 200 Hz |
| 每个输出端负载电容 | 5 pF |
| 电源电压范围 | VDD − VSS ≤ 40 V |

### 电极模型与 AC 耦合网络

| 项目 | 数值 / 要求 |
|---|---:|
| 皮肤-电极界面电阻 REL | 1 MΩ |
| 皮肤-电极界面电容 CEL | 10 nF |
| 电极模型 | REL–CEL tank |
| 输入读出方式 | AC-coupled |
| AC coupling network | 使用离散无源器件设计 |

### 工艺与器件约束

| 项目 |              数值 / 要求 |
|---|---------------------:|
| 设计环境 |     Cadence Virtuoso |
| TFT 工艺 |            AT_4000TG |
| 可用晶体管 |             pmos_TFT |
| TFT 沟道长度 L |   10 µm ≤ L ≤ 800 µm |
| TFT 沟道宽度 W | 50 µm ≤ W ≤ 50000 µm |
| 可用电容 |  MIM capacitor，`cap` |
| MIM 电容范围 |  500 fF ≤ C ≤ 200 pF |

### 鲁棒性与良率评估

设计需要在 typical、slow 和 fast 工艺角下进行验证，并结合 mismatch 进行 500-point Monte Carlo 仿真。

如果某个仿真点中一个或多个性能指标相对于目标规格偏差超过 20%，则该点被判定为 fail；否则判定为 pass。pass 点数占总仿真点数的比例定义为 yield。