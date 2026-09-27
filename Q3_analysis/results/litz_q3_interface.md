# litz_q3（队友 COMSOL 问题三）与 MATLAB PEEC 的对接评估

日期：2026-09-26 ｜ 来源：`Question/litz_q3`（仿真）、`Question/litz_topology_explanation`（几何验证）

---

## 1. litz_q3 到底是什么

**方法**：准三维（quasi-3D）= 连续三维中心线 + **COMSOL 二维截面 FEM 阻抗矩阵**（400×400）+ 逐轴向站位的股号置换电路，全长 1 m 串联求解。**不是完整三维体网格涡流。**

- 截面方程：$-\nabla\cdot(\nabla A_z/\mu_0)+j\omega\sigma A_z=\sigma E_j$，$J_z=\sigma(E_j-j\omega A_z)$；空气 $\sigma=0$，外边界 $A_z=0$。
- 损耗 $=\int |J|^2/\sigma\,dS$，**RMS 相量、不乘 1/2** —— 与我们的 PEEC/COMSOL 约定一致 ✅。
- 电路：$\bar Z\approx\langle P(z)^\top Z_{cs}P(z)\rangle$，$I=I_{tot}\bar Z^{-1}\mathbf1/(\mathbf1^\top\bar Z^{-1}\mathbf1)$；full_exchange 由循环对称性得 $I_j=20/N$。

**设计点（与我们的 PEEC 不同）**：$A_{Cu}=6.0$ mm²（我们 5.0953），$N$=400（主）/484（最终族）（我们 331），$d$=138.2 / 125.6 µm（我们 140），方形网格 20×20 / 22×22（我们六角环）；包络外径 6.587 mm；约束：间隙 ≥10 µm、编织角 ≤30°、弯曲半径 ≥5d、平均线长 ≤1.08。

**拓扑族**：`straight` / `rigid_twist`（刚性旋转，半径不变）/ `shell_exchange`（层内换位）/ `full_exchange`（偶数方格闭合 Hamilton 回路，逐股遍历全部站位）。

## 2. 主要结果（真实 COMSOL）

| 方案 | N | P / mm | Rac mΩ/m | Rdc mΩ/m | K=Rac/Rdc | P_cu W/m | η |
|---|---:|---:|---:|---:|---:|---:|---:|
| 初始整体扭转 | 400 | 100 | 13.806 | 2.895 | 4.769 | 5.523 | 901% |
| 推荐完整换位 | 484 | 500 | 4.398 | 2.938 | 1.497 | 1.759 | ≈0 |

相对初始：K −68.6%、铜耗 −68.1%。三种拓扑受控对比（N=400、P=250）：刚性 4.760 → 层内换位 4.663 → 完整换位 1.596（η 909%→358%→0）。验证：网格 −2.6e-6、外边界 −1.9e-8、单丝 Bessel 1.26e-5、矩阵互易 6.8e-16、矩阵预测 vs 直接 COMSOL 功率差 1.25e-12。

## 3. 数据接口清单

| 文件 | 内容 | 形状 |
|---|---|---|
| `data/final_484_centerlines.npz` | `xyz_m`；`strand_ids`；`topology`；`period_m` | (4001,484,3) [站位,股,xyz]，SI |
| `data/optimized_centerlines.npz` | 同上（N=400） | (1001,400,3) |
| `data/impedance_matrix.npz` | `impedance`(400×400 复)、`admittance`、`coordinates`(400,2)、`radius_m` | COMSOL 截面 |
| `data/recommended_circuit.npz` | `currents`(400,复)、`length_factors`(400,)、`impedance` | 阵后处理 |
| `data/raw/<case>/result.npz` | `currents`、`metrics`(6 行)、`strand_losses_W_per_m` | 每案例 |
| `data/*.csv/json` | topology_control、strand_count_control、fem_convergence、optimization_history(169 次)、final_results.json | 记录 |

`metrics` 6 行：数值铜面积 / 总复电流 / 截面 Rac / 截面 Rdc / 截面 Rac/Rdc / 20 A 截面铜耗。

## 4. 与 MATLAB PEEC 的对接：可行性与缺口

**兼容点** ✅
- 同一物理约定（RMS、$|J|^2/\sigma$ 无 1/2）；
- 两者都在"截面阻抗矩阵 + 置换/分段 + 并联共压"框架下，公式同族；
- 我们 `q3_solver.m` 的共压求解与他们的电路式数学等价。

**缺口** ⚠️
1. **npz → MATLAB**：MATLAB 不能原生读 `.npz`（非 HDF5）。WSL 侧 `python3` 有 numpy 2.4.6（无 scipy/h5py）→ 写个 `npz→CSV` 转换器即可。
2. **几何模型不同**：他们是"方格 Hamilton 回路 + 连续三次样条"，我们是"六角环 + 离散分位数层路径"。`q3_solver` 目前**自生成几何**，不接受外部中心线；需加一个适配器：外部中心线 → 按 $K$ 站位采样出 $(X,Y,V_X,V_Y,S)$ → 复用同一互感矩阵与共压求解。
3. **口径不同（关键）**：他们的 $K=R_{ac}/R_{dc,geom}$ 来自**含股内邻近效应的截面 FEM**（他们 straight N=400 就是 4.76）；我们的 $J$ 只含"股间+单丝趋肤"，$J_{full}$ 才补股内项。**对接口径：我们的 $J_{full}$ 对他们的 $K$；单独 $J$ 会系统性偏低。** 这也解释了为何我们模型里换位后 $J\to1.0$ 而他们的完整换位停留在 **1.497**——差额就是股内邻近损耗，只能由 COMSOL 级模型给出。
4. **设计点不同**：$A_{Cu}$ 6.0 vs 5.0953、$N$ 400/484 vs 331、$d$ 125.6/138.2 vs 140 µm。绝对 Rac 不可直接比；要对照必须锁定同一 $(N,A_{Cu},d)$。

## 5. 建议的对接动作（按优先级）

1. **转换器**：`Question/Q3_PEEC/litz_npz_to_csv.py`（numpy only）——把 `final_484_centerlines.npz`、`impedance_matrix.npz`、`recommended_circuit.npz` 导出为 CSV。
2. **MATLAB 适配器**：`q3_ingest_litz.m`——读中心线 CSV，按 $K$ 站位采样，**用同一 PEEC 核**复算逐股电流与 $R_{ac}$，与 `recommended_circuit.npz` 的 COMSOL 电流/`final_results.json` 的 4.398 mΩ/m 对照。
3. **同口径交叉验证**：在**同一设计点**（建议取他们 N=400 的三个拓扑）上，PEEC 与他们的截面 FEM 比 $R_{ac}$ 排序与逐股电流；确认我们的 $J_{full}$ 能贴到他们的 $K$。
4. **收口**：论文里把两种方法并列为"PEEC 快速扫描 + COMSOL 截面精算"，统一按 $K=R_{ac}/R_{dc,geom}$ 报告，并声明各自的建模近似。

## 6. 其他

- 他们完整换位 **1.497** 是对我们前面"完美换位地板"的**实测修正**（我们旧估 $J_{full}\sim$1.7–2.4 偏保守）；说明把 $R_{ac}/R_{dc}$ 压到 1.5 附近是可信目标，继续压"股间通道"无空间，**只能靠减 $d$ / 改装填**（与 §8 结论一致）。
- `Question/第三题_Matlab.zip` 内含另一路 agent 的 `Q3_PEEC/q3_calibrate.m`；该脚本记录了我们清理旧核心（`q3_peec_core.m`）这一并发事件，并做了"存在则用旧核心、否则退用 `q3_solver.m`"的兜底。与 litz_q3 无关，属另一条线。

---

## 7. 临时转换与 PEEC↔COMSOL 对接结果（2026-09-26 实测）

**转换**：`Q3_PEEC/litz_npz_to_csv.py`（numpy only）导出中心线（~500 站位）与 COMSOL 逐股电流为 CSV。

**方法**：`Q3_PEEC/q3_ingest_litz.m` 读中心线 → 有限差分切线/线长因子 → 同一 PEEC 互感矩阵与共压求解（Bessel 内阻抗 + $\ln$ 互感），与 COMSOL 发布值对照。

| 方案 | N | $K_{peec}$ | $K_{comsol}$ | 偏差 | $\eta_{peec}$ | $\eta_{comsol}$ |
|---|---:|---:|---:|---:|---:|---:|
| rigid_twist (P=0.25) | 400 | 4.676 | 4.761 | **−1.77%** | 905% | 909% |
| shell_exchange (P=0.25) | 400 | 4.571 | 4.663 | **−1.97%** | 355% | 358% |
| full_exchange (P=0.25) | 400 | 1.117 | 1.596 | −30.0% | 1.6% | ≈0 |
| final_484 (P=0.5) | 484 | 1.040 | 1.497 | −30.5% | 1.3% | ≈0 |

逐股电流（final_484，COMSOL verified vs PEEC）：$|I|$ 最大偏差 **1.08%**、复数最大偏差 1.32%（相对 $\bar I$）、相关系数 **0.99999**。

**判读**：
1. **电流分布与未换位/弱换位损耗被高精度复现**（rigid/shell 偏差 ~2%，$\eta$ 与 COMSOL 几乎相同）——证明我们的互感 PEEC 与他们的截面 COMSOL 在同一几何上一致。
2. **完整换位方案 PEEC 系统性低 30%**，差额正是 PEEC 未解析的**股内邻近损耗**（均流后成为主导）。这与 §4 的口径判断完全吻合：**我们 $J$（股间+趋肤）↔ 他们的 $K$ 差一个股内通道**。用他们的 COMSOL 定量：full_exchange 股内邻近使 $K$ 从约 1.04 抬到 1.497（$P$ 从 1.22 抬到 1.76 W/m，约 +0.54 W/m）。
3. 因此最终报告口径应统一：**$R_{ac}/R_{dc}$ 的"完整值"以 COMSOL 截面为准（full_exchange≈1.50），PEEC 用于方案排序与股间通道快速评估**。

**改动记录**：临时转换器 `litz_npz_to_csv.py`、ingester `q3_ingest_litz.m`、结果 `litz_peec_vs_comsol.csv` 均在 `Q3_PEEC/`（中心线 CSV ~8–11 MB，属临时产物，可随时重生成）。

