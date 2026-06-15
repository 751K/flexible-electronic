# 技术文档

## 目录

- [电路拓扑](#电路拓扑)
- [器件模型](#器件模型)
- [DC 求解器](#dc-求解器)
- [AC 小信号求解器](#ac-小信号求解器)
- [噪声求解器](#噪声求解器)
- [叉指数（NF）](#叉指数nf)
- [工艺角与失配](#工艺角与失配)
- [校准方法](#校准方法)
- [Cadence Virtuoso 工程说明](#cadence-virtuoso-工程说明)

---

## 电路拓扑

AFE 是一个 **10 管全差分放大器**，包含两级：

```
第一级：差分输入对
  M6   — 尾电流源（VB 偏置）
  M7/8 — 输入差分对（VCM 输入，drain 接 VOP/VON）

第二级：交叉耦合正反馈电平移位器
  M11    — 电平移位器尾电流（VC 偏置）
  M12/13 — 交叉耦合正反馈对（gate 接 VOP/VON，drain 接 VFBN/VFBP）
  M14/15 — 电平移位器负载（gate = GND）

输出级：
  M9/10  — 共源放大（gate 接 VFBP/VFBN，drain 接 GND）
```

**MNA 求解节点**（6 个）：

| 索引 | 节点 | 说明 |
|------|------|------|
| 0 | VOP | 正输出 |
| 1 | VON | 负输出 |
| 2 | VFBP | 交叉耦合输出（正侧） |
| 3 | VFBN | 交叉耦合输出（负侧） |
| 4 | NET20 | 电平移位器尾节点 |
| 5 | NET2 | 输入对尾节点 |

**偏置电压**：VDD = 40V，VCM ≈ 30.65V，VB ≈ 9.84V，VC ≈ 16V。

拓扑定义在 `topology.py` 的 `AFE_TOPO` 中，是「单一真源」——DC KCL、逐管偏置映射、AC/噪声端子表全部从器件表派生，不在求解器里重复手写。

---

## 器件模型

`pmos_tft_model.py` 是 **AT_4000TG PMOS-OTFT** 的 Verilog-A 等价 Python 实现。

### 模型结构

器件内部包含 **接触网络节点** s1/d1，在 source/drain 端子与沟道之间引入接触电阻和接触肖特基模型：

```
Source (s) ──[contact]── s1 ──[internal R]── d1 ──[channel + leak]── Drain (d)
```

DC 求解时，先对内部节点 s1/d1 做 `fsolve` 收敛，再计算端子电流 `get_Idc`。

### 沟道模型

基于 **Vissenberg-Matters (VM)** 迁移率模型，采用 percolation theory 描述的 disordered organic semiconductor：

- 载流子迁移率服从指数温度依赖：μ ∝ (T/Tt)^(Tt/T)
- 电流以 `Vss * softplus((V-Vg+Vfb)/Vss)` 的形式表达，避免亚阈值区的不连续
- 包含 channel length modulation（λ）

### 噪声模型

- **Thermal noise**：`4kT·gm·2/3` + `2q·(Ich+Ioff)`
- **Flicker noise (1/f)**：Hooge 模型，`S_flicker = (Hooge·q·Ich²) / (|W·L·Ci·(Vd1-Vg+Vfb)|) / f`

### 寄生电容

Cgs/Cgd 由三部分叠加：边缘电容（overlap）、沟道电容（arcian 过渡）、面积电容。通过 `get_capacitances` 返回。

### 关键参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| VT | −3.03 V | 阈值电压 |
| Ci | 2.4 | 单位面积绝缘层电容 |
| NF | 1 | 叉指数 |
| Reg | 1 | 接触电阻模式 |
| alfa | 4.5455×10⁷ | 有效态密度 |
| Nt | 1×10²¹ | 局域态密度 |
| Tt | 305 K | 特征温度 |

---

## DC 求解器

### 基本方法

`ac_solver.py` 中的 `ac_solve` 对 6 节点 KCL 方程组用 `scipy.fsolve` 求解。KCL 残差由 `topology.py` 的 `dc_residuals` 方法从拓扑表自动派生：每个器件的 drain 贡献 +Id，source 贡献 −Id，再加 gmin·V 项。

### 鲁棒性三守护

多稳态电路（交叉耦合正反馈结构）需要确保收敛到 Spectre 对应的物理分支。求解器有三层守护，只在冷启动/异常时触发：

1. **源步进连续法**（Source-Ramp Continuation）—— 主种子策略
   - 在对称系统（4 未知数）上，将所有偏置电压从 0 线性 ramp 到目标值
   - 追踪上电过程的物理轨迹，模拟 Spectre 的 pseudo-transient 行为
   - 19 步 ramp（0.1→1.0），每步以 fsolve 跟踪解

2. **对称守护**（Symmetry Guard）
   - 无失配时物理解应对称（VOP=VON, VFBP=VFBN）
   - 若 fsolve 锁存到对称破缺的假根，用对称平均作为种子重解 4 节点对称系统

3. **物理性守护**（Physicality Guard）
   - 节点电压超出 [0, VDD] 即为非物理分支
   - 以轨内种子重解，选轨内解

### 快速路径

带种子（`x0_guess` 参数）调用时跳过连续法，直接走 fsolve。MC/优化等 in-loop 扫描利用名义解作种子，速度大幅提升。

---

## AC 小信号求解器

### MNA 构建

每个频率点，构建 6×6 复 Y 矩阵，对每个器件做 `_stamp_mos`（`ac_mna.py`）：
- **跨导**（gm）：VCCS stamp，id = gm·(Vg−Vs)，注入 drain/source
- **输出电导**（gds）：drain-source 间导纳
- **电容**（Cgs, Cgd）：gate-source/gate-drain 间 sC 导纳

输入信号作为已知电压源出现在 RHS：`M7.gate = +0.5V, M8.gate = −0.5V`（差模 1V）。

5pF 负载电容加在 VOP/VON 对 GND。

### 端子 gm/gds

**关键差异**：AC/噪声用的 gm/gds 是 **端子值**，通过对 `get_Idc` 做有限差分（ΔV = ±1mV）提取，而非器件内部的沟道 gm。

原因：OTFT 的接触电阻网络会 degenerate 沟道 gm。用沟道 gm 时增益偏差 0.8 dB、BW 偏差 18 Hz；用端子 gm 后偏差 <0.05 dB、<0.1 Hz，与 Spectre 一致。

### 输出

- 增益曲线（|VOP−VON| / vin_diff 扫频）
- DC 增益（dB）
- −3dB 带宽
- 全节点 DC 工作点
- 每管小信号参数（gm, gds, Cgs, Cgd）

---

## 噪声求解器

### 方法（与 Spectre 一致的转移阻抗法）

`noise_solver.py` 在与 AC 相同的 6 节点 MNA 上做噪声传播：

1. **DC 求解**：复用 `ac_solve` 得工作点和小信号参数
2. **构建 Y 矩阵**：与 AC 完全一致（输入对 gate 为 AC ground，drive={}）
3. **每管噪声注入**：对每个晶体管，在其 drain 和 source 之间注入单位电流 I_inj
   - 求解 Y·V = I_inj，读差分输出电压 (VOP−VON)
   - 转移阻抗 Z_k = (VOP−VON) / I_inj
   - 该管的输出噪声贡献 PSD = |Z_k|² × S_id,k
4. **总和**：各管噪声不相关 → 直接叠加 PSD → 总输出噪声 PSD
5. **折合输入**：除以放大器增益 |H(f)|² → 输入参考噪声 IRN PSD

### 噪声注入位置

注入在器件的 drain/source 节点（从 AC 端子表派生），对于：
- **rail 端**（如 M9 drain = GND）：AC ground → 不贡献 Z（Z = 0）
- **solved 端**：Z = Yinv[VOP, node] − Yinv[VON, node]

### 积分

`band_rms` 用 `scipy.trapezoid` 对 PSD 在指定频带内积分再开方，得 RMS 电压。

---

## 叉指数（NF）

### 作用范围

NF（Number of Fingers，叉指数）影响器件的**栅寄生电容** Cgs/Cgd，因为叉指几何决定了 fringe/overlap 电容的面积和边缘周长。NF **不影响** DC 电流或 gm，因为电流由总 W/L 决定。

具体影响链：

```
NF ↑ → 叉指数增加 → fringe/overlap 电容变大 → Cgs/Cgd ↑
                                                   ↓
                                    AC: −3dB BW 降低（电容负载增大）
                                Noise: 高频段电容耦合增加
```

### 使用方式

与 `corner` 相同的三种传入模式：

```python
# 不传（默认 NF=1）
r = noise_analysis(sizes, bias, freqs)

# 全局 int（所有管相同叉指数）
r = noise_analysis(sizes, bias, freqs, nf=4)

# 逐管 dict（每个管独立 NF，缺省 = 1）
nf_map = {"M7": 120, "M8": 120, "M9": 8, "M10": 8}
r = noise_analysis(sizes, bias, freqs, nf=nf_map)
```

### 实现链路

`_dev_nf(nf, name)` 解析 NF → 传入 `PMOS_TFT(NF=nf)` → 模型内部 `_precompute_constants` 根据 NF 重算 `fw = W/NF` 和栅电容几何参数 → `get_capacitances` 返回修正后的 Cgs/Cgd → AC 求解器和噪声求解器自然使用修正后的电容值。

---

## 工艺角与失配

求解器支持两种不理想因素，通过 `corner` 参数统一透传：

### 全局工艺角（Process Corner）

```python
slow  = {"pvt0": -3*0.0753, "pbeta0": -15*0.036}
fast  = {"pvt0": +3*0.0753, "pbeta0": +15*0.036}
```

所有管施加相同偏移。

### 逐管失配（Mismatch）

```python
mc = {"M7": {"mvt0": σ_vt_7, "mbeta0": σ_β_7}, ...}
```

每管独立偏移量。σ 为：
- σ_vt = `Avt / sqrt(W·L)`，Avt ≈ 0.0753 V·µm
- σ_β = `Aβ / sqrt(W·L)`，Aβ ≈ 0.036

### 参数作用

| 参数 | 作用 | 模型内部影响 |
|------|------|-------------|
| pvt0 | 全局 VT 偏移 | Vfb = VT·(Ci_ratio)·(1+pvt0) |
| mvt0 | 局部 VT 偏移 | mvt0/sqrt(W·L) 附加到 Vfb |
| pbeta0 | 全局 β 偏移 | beta 缩放 (1+pbeta0) |
| mbeta0 | 局部 β 偏移 | beta 缩放 (1+mbeta0) |

---

## 校准方法

求解器针对 **Cadence Spectre 24.1.0.078** 校准，方法如下：

### 校准步骤

1. **DC 工作点**：对照 Spectre DC 仿真，确保所有 6 节点电压一致（含内部节点）
2. **AC 增益/带宽**：对照 Spectre AC 扫频，确认增益曲线重合
3. **噪声 PSD**：对照 Spectre noise 分析，确认输出噪声谱和每管贡献比例
4. **工艺角**：对照 Spectre 工艺角仿真（typical/slow/fast）
5. **Monte Carlo**：对照 Spectre MC 仿真（500 次逐管失配），比较均值与 σ

### 关键校准点

- **端子 gm** 而非沟道 gm（消除接触电阻退化的系统性偏差）
- **接触模型** 极性排序（`_va_sorted_nodes` 模拟 Verilog-A 的 ternary operator）
- **Hooge 系数** 0.05（噪声幅度匹配 Spectre）
- **gmin** 1e−12（同 Spectre 默认 gmin 值）

### 校准精度

| 项目 | 偏差 |
|------|------|
| 典型 gain | ±0.01 dB |
| 典型 BW | 同点 |
| 典型 IRN | ≤ 数% |
| 工艺角 gain/BW/IRN | 同上 |
| MC 均值 | 一致 |
| MC σ | 一致 |
| 全范围随机点（含多稳态/退化/超轨边缘）| gain 0.01 dB / IRN 2% |

---

## Cadence Virtuoso 工程说明

### 单元库

| 目录 | 内容 | 说明 |
|------|------|------|
| `FGB/` | schematic + layout + symbol + av_extracted | 完整单元库，含提取后版图 |
| `FGB_layout/` | schematic + layout + symbol | 版图版本单元库 |
| `Two/` | schematic + symbol | 两管结构辅助单元 |
| `top/` | schematic + symbol | 顶层原理图 |
| `top_tb/` | schematic + maestro | 测试台 + Maestro 仿真配置 |

### 打开方式

在 Cadence Virtuoso 中，将项目根目录添加到 Library Path，即可浏览原理图、版图和仿真结果。

`top_tb/maestro/` 下保存了 Maestro 仿真配置和历史数据。
