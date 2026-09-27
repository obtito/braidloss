# Q3 优化变量选择：六判据的定量依据（纯标准库，可复现）
# 运行: python q3_variable_selection.py
# 数据来源（出处均在注释标明）:
#   Question/Q2_Circular_COMSOL/results_summary.json + ring_statistics.csv
#   Question/Q2_RegularTwist_COMSOL/results_summary.json
#   Question/PEEC/peec_vs_comsol_summary.csv
import math

# ---------- 统一工况与 Q2 基准 ----------
mu0 = 4e-7 * math.pi
sigma = 5.8e7
f = 200e3
w = 2 * math.pi * f
I = 20.0
delta = math.sqrt(2 / (w * mu0 * sigma))
d = 0.14e-3
a = d / 2
N = 331
Acu = N * math.pi * a ** 2
Pdc = I * I / (sigma * Acu)
p = 0.15e-3
R_env = 1.575e-3
fill = Acu / (math.pi * R_env ** 2)   # 0.6538

print("=" * 68)
print("Q3 VARIABLE SELECTION - quantitative basis for each criterion")
print("=" * 68)
print(f"delta = {delta*1e6:.2f} um ; A_Cu = {Acu*1e6:.4f} mm^2 ; P_dc = {Pdc:.4f} W/m")

# ---------- [1] 度量护栏：J 对丝长因子 s_bar 一阶不变 ----------
print("\n[1] Metric guarantee: J is first-order invariant to strand lengthening")
Rac_u, Rdc_u = 14.899678, 3.383748       # 圆形未换位 (Q2_Circular)
Rac_t, Rdc_t = 14.401182, 3.434493       # 规则绞合 P=40mm (Q2_RegularTwist)
s_bar = Rdc_t / Rdc_u
J_u, J_t = Rac_u / Rdc_u, Rac_t / Rdc_t
Rac_cf = Rac_u * s_bar                    # 反事实：纯几何拉长、无电磁改变
print(f"  s_bar = {s_bar:.5f}")
print(f"  counterfactual (geometry-only): J = {Rac_cf/Rdc_t:.4f} == J_untwisted {J_u:.4f} (构造性相等)")
print(f"  => 纯拉长不改变 J；绞合的 J 变化 {(1-J_t/J_u)*100:.2f}% 全部来自电磁/分布效应")
print(f"  但每米 Rac {'下降' if Rac_t<Rac_u else '上升'} {abs(1-Rac_t/Rac_u)*100:.2f}% 且铜重/m 上升 {(s_bar-1)*100:.2f}%")
print("  => J 无法惩罚'过度倾斜'；必须同时报告 P(W/m) 与铜重/米（防度量博弈）")

# ---------- [2] 换位的代价：股内通道上升（变量分层依据） ----------
print("\n[2] Transposition trade-off: intra-strand channel RISES (variable tiering)")
C = sigma * math.pi * a ** 4 * w ** 2 / 8  # W/m per strand per T_rms^2
ring_I = [0.00054018, 0.00056376, 0.00080732, 0.00148707, 0.00293551,
          0.00596772, 0.01239393, 0.02612674, 0.05566997, 0.11959078,
          0.25856026]                      # COMSOL ring_statistics.csv 逐环均值
cnt = [1, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60]
B2w = 0.0
for k in range(11):
    if k == 0:
        B = 0.0
    else:
        r = k * p
        Iin = sum(cnt[j] * ring_I[j] for j in range(k))
        Iown = cnt[k] * ring_I[k]
        B = mu0 / (2 * math.pi * r) * (Iin + Iown / 2)
    B2w += cnt[k] * B * B
P_intra_un = C * B2w
B2_eq = (mu0 * I) ** 2 / (8 * math.pi ** 2 * R_env ** 2)
P_intra_eq = C * N * B2_eq
print(f"  untwisted (imbalanced) : <B^2> = {B2w/N*1e6:.2f} mT^2 -> P_intra = {P_intra_un:.3f} W/m (COMSOL 0.396)")
print(f"  perfect transposition  : <B^2> = {B2_eq*1e6:.2f} mT^2 -> P_intra = {P_intra_eq:.3f} W/m")
print(f"  ratio x{B2_eq/(B2w/N):.2f} -- 换位唤醒内部闲置股, 平均磁场暴露反而上升")
print(f"  net: harvest 4.211 - pay {P_intra_eq-P_intra_un:.2f} = {4.211-(P_intra_eq-P_intra_un):.2f} W/m")
print(f"  floor J = 1 + {P_intra_eq/Pdc:.3f} = {1+P_intra_eq/Pdc:.3f}")
print("  => 变量按'作用区间'分层: 换位前股间主导(S,alpha); 换位后地板主导(d,排布)")
d_alt = 0.10e-3                          # 敏感性轴: d 0.14 -> 0.10, 保持 A_Cu
N_alt = Acu / (math.pi * (d_alt / 2) ** 2)
P_intra_alt = 0.396 * (d_alt / d) ** 2   # 固定 A_Cu: 总股内损耗 ∝ N·a^4 ∝ d^2
floor_d14 = 1 + fill * N * (d / delta) ** 4 / 256
floor_d10 = 1 + fill * N_alt * (d_alt / delta) ** 4 / 256
print(f"  d 效应两副面孔: 未换位基线上仅省 {(0.396 - P_intra_alt) / 5.960 * 100:.1f}% "
      f"(~可辨识阈值 2.5%, 排不动) vs 换位后地板 {floor_d14:.3f} -> {floor_d10:.3f} "
      f"(超额 0.681 -> {floor_d10 - 1:.3f}, 近乎减半, 远超噪声)")

# ---------- [3] alpha 的两翼机制 ----------
print("\n[3] Braid angle alpha - two-sided mechanisms:")
R_mm = R_env * 1e3
print(f"  {'a/deg':>6} {'cos2a':>7} {'1/cos a':>8} {'cross/m(8载)':>13}")
for adeg in [13.3, 20, 25, 30, 35, 40, 45]:
    ar = math.radians(adeg)
    cross = 16 * 1000 * math.tan(ar) / (math.pi * R_mm)  # 16 反向族对, 每对/m
    print(f"  {adeg:6.1f} {math.cos(2*ar):7.3f} {1/math.cos(ar):8.3f} {cross:13.0f}")
print("  cos2a   = 交叉点处反向族互感的切线余弦 (45 度时局部互感消失; 平面电感解耦原理)")
print("  1/cos a = 铜重/每米 (J 一阶豁免, 但质量与 W/m 惩罚)")
print("  cross/m = 交叉点密度 16*tan(a)/(pi*R)*1000 (管状编织 8 载纱器/16 反向对)")
a_q2 = math.degrees(math.atan(2 * math.pi * 1.5e-3 / 40e-3))
print(f"  Q2 同向绞合外层 a = {a_q2:.1f} deg (tan a = 2*pi*r/P = 0.2356)")
print(f"  => Q2 实测 EM 收益 4.77% 即小 alpha 边缘的 O(a^2) 修正下限")

# ---------- [4] 子束粒度 m：ln 核精确计算 ----------
print("\n[4] Sub-bundle granularity m: internal imbalance from ln-kernel (hex)")

def hex_pos(nr):
    pos = [(0.0, 0.0)]
    for k in range(1, nr + 1):
        for j in range(6 * k):
            th = 2 * math.pi * j / (6 * k)
            pos.append((k * p * math.cos(th), k * p * math.sin(th)))
    return pos

def dX_ln(pos, Itot):
    n = len(pos)
    Is = Itot / n
    def sum_ln(px, py, excl):
        s = 0.0
        for i, (x, y) in enumerate(pos):
            if i == excl:
                continue
            s += math.log(math.hypot(x - px, y - py))
        return s
    kmax = (n - 1) // 6
    ei = 1 + 6 * (kmax - 1)          # 最外环第一股
    ex, ey = pos[ei]
    dL = (mu0 / (2 * math.pi)) * Is * (sum_ln(ex, ey, ei) - sum_ln(0, 0, 0))
    return w * dL

Rs = 1 / (sigma * math.pi * a ** 2)
print(f"  {'m':>4} {'I_sub/A':>8} {'dX_sub':>9} {'dX/Rs':>8} {'eta上限':>8}")
for nr in [1, 2, 3, 4, 5]:
    pos = hex_pos(nr)
    m = len(pos)
    Isub = I * m / N
    dX = dX_ln(pos, Isub)
    print(f"  {m:4d} {Isub:8.3f} {dX:9.4f} {dX/Rs:8.4f} {dX/(Rs+dX)*100:7.2f}%")
print("  => 子束越粗, 束内残余失衡越大: m=19(8%) ~ m=37(16%) ~ m=61(25%)")
print("     与束级 dX=2.51 vs R_sub=0.042 的尺度论证一致: 每一级都要换位")

# ---------- [5] 可辨识性门槛 ----------
print("\n[5] Identifiability thresholds (from Q2 PEEC-vs-COMSOL statistics):")
print("  |PEEC-COMSOL| gap: hex 3.9%, circular 3.65%, twist 2.54% (absolute)")
print("  Delta reproduction error (twist vs circular): 0.177 mOhm = 1.2% of Rac")
print("  => PEEC 可排序门槛: 效应 >= 2.5%; d/排布类小效应必须 COMSOL 定案")
rows = [
    ("S 换位方案族", "股间 71%", "4.40 -> 1.7~2.5 (-50%)", "主变量 (PEEC)"),
    ("alpha 编织角", "股间+耦合", "两翼机制, 45deg 处 O(1) 解耦", "主变量 (PEEC)"),
    ("m 子束粒度", "股间残余", "m=19~61: 8%~25% 失衡上限", "主变量 (PEEC)"),
    ("d 丝径", "股内地板 6.6%", "基线上 3.3% / 换位后 1.68->1.35", "敏感性 (COMSOL)"),
    ("排布 盘/环", "股内地板", "地板 1.68 -> 1.14 (R=3mm)", "第二层 (COMSOL)"),
    ("N 股数", "与 d 锁死", f"Nd^2 = {N*d**2*1e6:.4f} mm^2 常数", "协议固定"),
    ("p 节距", "封装约束", "~0", "协议固定"),
]
print(f"  {'变量':<14}{'主通道':<12}{'对 J 的效应量':<30}{'裁决'}")
for r in rows:
    print(f"  {r[0]:<14}{r[1]:<12}{r[2]:<30}{r[3]}")

# ---------- [6] 多级结构整除性 ----------
print("\n[6] Multi-level feasibility (N=331 is prime!):")
print(f"  {'N':>5} {'d/mm':>7} {'dd%':>6} {'分解':<10} {'m<=7自平衡':<12} {'m in 8..32需绞'}")
for NN in [331, 320, 336, 352, 385]:
    dd = math.sqrt(4 * Acu / (math.pi * NN))
    fac = []
    x = NN
    for q in range(2, 20):
        while x % q == 0:
            fac.append(q)
            x //= q
    small = [m for m in range(5, 8) if NN % m == 0]
    mid = [m for m in range(8, 33) if NN % m == 0]
    print(f"  {NN:5d} {dd*1e3:7.4f} {(dd/d-1)*100:+6.1f} {'x'.join(map(str, fac)) or '(素数)':<10} "
          f"{str(small):<12} {mid}")
print("  => 多级/编织优先可整除 N: 320(d+1.7%) / 336(d-0.7%), 且 336 = 48x7 恰好给出 m=7 自平衡子束")
print("     352/385 的 d 偏离 3~7%: 换位地板按 d^4 进一步改善, 但漆膜占比与工艺需重新校验, 需单独说明")

# ---------- [7] 最终参数化 ----------
print("\n[7] Final parameterization:")
print("  主变量   : S(方案族), alpha(编织角/周期), m(子束粒度)")
print("  协议固定 : d=0.14, N=331, p=0.15, A_Cu, f, I, L (受控比较)")
print("  敏感性轴 : d in {0.10, 0.14}; 排布 in {实心盘, 环 R=3mm}")
print("  报告协议 : J, P(W/m), 铜重/米, 包络 D, eta_I, 停留分布偏差")
