# 问题四：最优拓扑逼近方案与 MATLAB 脚本设计

日期：2026-09-26 ｜ 依据：`Question/Q4_analysis/问题四_完美换位判据与拓扑设计.md`（v2）§4 设计族、§5 交接规范（T1–T5）。

---

## 0. 从判据到逼近方案的一步关键认识

v2 判据：**边际公平（主）** $O_{im}=1/N$ 可精确达到；**股对公平（加强）** 一般只能逼近；外层再最小化 $\Phi=\sum_i\int|B_i|^2dz$。

**关键认识（本方案基础）**：v2 §2.3 要求"等面积站位"。对六角环等间距布局（第 $k$ 环半径 $r_k=kp$、$6k$ 个角位），角向间距 $2\pi r_k/(6k)=\pi p/3$ **与 $k$ 无关**，故每站 Voronoi 面积 $\approx(\pi p/3)\cdot p$ 近似相等；环 $k$ 面积 $\propto 6k$。因此 **$p_k=6k/N$ 的壳容量就是"等面积站位"的离散实现**，与 §2.3 路径一相容——`q3_solver.m` 的分位数层路径天然满足该约定。这把"等面积"从额外的几何构造简化成了一个已有性质。

由此，逼近方案 = 在 `q3_solver` 已有的**精确公平**基础上，做两件事：
1. **证书化**：算 $D_O$（边际公平）、$D_E$（股对公平）作为"完美度"证据；
2. **压 $\Phi$**：在保持公平的前提下改变**装填**（空心化 / 密度渐变 / 减 $d$）与**周期/角**（$P,\alpha$）使残余损耗最小。

---

## 1. 逼近方案（设计族 $\mathcal F$ 的具体化）

| 维度 | 选项 | 取值/说明 |
|---|---|---|
| 站位集 | 六角环等间距（≈等面积，见 §0） | 环 $k$：$6k$ 站，$r_k=kp$；$N=1+\sum 6k$ |
| 半径路径 $\rho$ | 分位数三角波（等时=等面积） | $r_i(z)=F^{-1}(\mathrm{tri}(\varphi_i+z/P))$，$D_O=0$ |
| 角向调度 | 刚体旋转 $\beta=\tan\alpha/r_{\max}$ + 层间槽位重排 | 占位守恒 ⇒ 间距 $\ge2a$ 由构造保证 |
| 周期 | $P\in[20,240]$ mm | 扫描（问题三：单调降并饱和） |
| 编织角 | $\alpha\in[2°,45°]$ | 扫描（问题三：单调升，取小角） |
| **装填变体** | 实心 / **空心环** / 密度渐变 | 新维度，Q4 的核心增量 |
| 丝径 | $d$ 敏感性 | 等 $A_{Cu}$ 下随 $N$ 联动（$\beta:0.95\to0.6$） |

**装填变体的实现（保持 $N,A_{Cu},$ 最小间距不变）**：
- 空心环：删去中心 $k_0$ 层，把其股数按等面积重分布到外层，外层环数增加、环间距略缩（仍 $\ge2a$）；内腔面积占比 $\eta_{\rm hollow}>0$；
- 密度渐变：环站数权重改为 $w_k\propto r_k^{1+\gamma}$，$\gamma>0$ 外密内疏。

**$\Phi$ 代理（不依赖场求解，用于排序）**：
$$
\Phi_{\rm proxy}=\underbrace{\text{装填均匀度偏差}}_{\text{站面积方差}}+\lambda_1\,\underbrace{\text{内腔面积占比}}_{\text{空心奖励}}-\lambda_2\,\underbrace{\max|dr/dz|}_{\text{曲率惩罚}} .
$$
Top 候选再用 PEEC 复算真实 $J$（`q3_solver` 核）；最终 $K$ 交 COMSOL。

---

## 2. MATLAB 脚本设计（对齐 §5 交接规范）

### 2.1 文件与职责

| 文件 | 职责 | 产物 |
|---|---|---|
| `Q4_PEEC/q4_peec_core.m` | 从 `q3_solver.m` 抽出的**通用核**：输入任意 $(X,Y,V_X,V_Y,S,a)$ → 共压求解 → $J,R_{ac},P,\eta_I,\tilde I_i$ | 复用函数 |
| `Q4_PEEC/q4_topology.m` | 生成器：cfg → 轨迹族、$s_i$、事件表、$D_O$、可行性 | **T1** `q4_<name>_centerlines.csv`、**T2** `q4_<name>_events.csv`、**T4** 可行性行 |
| `Q4_PEEC/q4_certificates.m` | 汇总完美度证书：$D_O$（CSV 统计）、$\eta_I/J$（PEEC）、$D_E$（有场则算，否则标注待补） | **T3** `q4_certificates.csv` |
| `Q4_PEEC/q4_run.m` | 枚举设计族、$\Phi$ 代理排序、Top3 交 PEEC、Top1–2 出 COMSOL 交接卡、双通道日志 | **T5** 排序表 + `q4_optimal.json` |
| `Q4_PEEC/q4_certify_existing.m` | 对**队友已有中心线**（484 full_exchange）与我们的 braid 直接算 $D_O$，验证判据 | T3 对既有方案的证书 |

### 2.2 接口与复用（避免重复造轮子）
- `q3_solver.m` 的共压求解段**原样抽成** `q4_peec_core`，`q3_solver` 改为调用它，保证 Q3/Q4 数值一致；
- 外部中心线路径沿用 `q3_ingest_litz.m` 的读入/切线/线长逻辑（已与 COMSOL 验证）；
- 几何生成沿用分位数层路径 + 槽位重排；装填变体只在**布局生成**处改，不下沉到求解核。

### 2.3 输出契约（逐项对应 v2 §5）
- **T1**：中心线 CSV，列 `strand,z,r,theta,s,section`，采样 $\le P/64$；
- **T2**：Hamilton/槽位事件表，列 `section,strand_from,strand_to,sigma_k,parallel_group`；braid word 展开 $\Delta=\sigma_1\cdots\sigma_{N-1}$、一周期 $=\Delta^{N}$；
- **T3**：`q4_certificates.csv`，列 `topology,N,d_mm,P_mm,alpha_deg,packing,D_O,D_E,eta_I,J_peec,K_comsol,note`；
- **T4**：`min_center_dist_um, min_bend_diameters, sbar_max, feasible, rejection`；
- **T5**：$\Phi$ 代理排序 + Top3 PEEC $J$ + Top1–2 交接卡（`scheme_card.json` + `centerlines.csv`）。

### 2.4 阶段门
| 门 | 判据 |
|---|---|
| G0 | `q4_peec_core` 在 straight/twist 上逐位复现 `q3_solver`（14.357230 / 14.035939） |
| G1 | 生成的公平拓扑 $D_O<10^{-3}$（构造精确，回归断言） |
| G2 | 可行性：$\min$ 间距 $\ge2a$、弯曲半径 $\ge5d$、$\bar s<1.02$ |
| G3 | 排序：$\Phi$ 代理排序与 PEEC $J$ 排序一致（Kendall $\tau$ 检查） |
| G4 | ≥1 拓扑满足验收（$D_O<10^{-3}$、$\eta_I<2\%$、可行），且相对 Q3 最优有量化 $K$ 下降预期 |

---

## 3. 预期与诚实边界

- **公平性**：$D_O\approx0$ 由构造保证（不是优化出来的），预期 $10^{-3}$ 以下；这是"完美换位"的硬证据。
- **残余 $K$**：地板由装填决定。参考 teammate 实测 $K=1.497$（N=484 实心方格）与预研地板公式（$d$、空心环）；本方案用**减 $d$ 与空心化**在装填层压 $\Phi$，给出 $K$ 下降的机理与预估。
- **口径警告（承接 workbuddy 记录）**：股内邻近项常量与均匀电流场存在口径问题（$/4$ vs $/8$、均匀电流低估），$D_E$ 与 $K$ 的**绝对值**须用重新推导后的口径；**排序与 $D_O$ 不受影响**。
- 全部结论需 PEEC/COMSOL 数据支撑；$\Phi_{\rm proxy}$ 只用于筛选，不作最终数值。

---

## 4. 待确认决策

1. **设计点 $N$**：(a) 331（延续 Q3，六角 $m$=10）；(b) 469（六角 $m$=12，≈队友 484）；(c) 直接用队友 484 方格布局。建议 (a)+(b) 两档，便于与 Q3 和 teammate 双向对照。
2. **$D_O$ 口径**：站（slot）等面积（推荐，与分位数构造一致）还是壳层面积权重；两者都算并报告。
3. **离散 braid word 版**：是否要在连续分位数版之外，额外实现"偶格 Hamilton + 样条"离散版（T2 事件表更完整、更贴制造语言）；建议连续版先行，离散版作为 T2 增强。
