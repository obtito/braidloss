# Q3 预研：机理分析与设计变量选择（全部数值可复现，纯标准库）
# 运行: python q3_preliminary_analysis.py
# 依赖: 无（刻意只用 stdlib，便于在任何环境复现）
import math

# ---------------- 统一工况 ----------------
mu0 = 4e-7 * math.pi
sigma = 5.8e7          # S/m
f = 200e3              # Hz
w = 2 * math.pi * f
I = 20.0               # A rms
delta = math.sqrt(2 / (w * mu0 * sigma))

# ---------------- Q2 基准设计 ----------------
d = 0.14e-3            # m
a = d / 2
N = 331
Acu = N * math.pi * a ** 2          # 5.0953 mm^2
Pdc = I * I / (sigma * Acu)         # 1.3535 W/m
R_strand = 1 / (sigma * math.pi * a ** 2)   # 1.12 ohm/m
p = 0.15e-3            # strand pitch
R_env = 1.575e-3       # bundle envelope radius (10 rings)
fill = Acu / (math.pi * R_env ** 2)  # 0.654

print("=" * 64)
print("Q3 PRELIMINARY ANALYSIS - mechanism & design variables")
print("=" * 64)
print(f"delta = {delta*1e6:.2f} um ; A_Cu = {Acu*1e6:.4f} mm^2 ; fill = {fill:.4f}")

# ---------------- [1] 两个失效模式独立量化 ----------------
# 单丝趋肤: Ferreira low-freq expansion, R_ac/R_dc ~ 1 + (1/48)(a/delta)^4
skin_excess = (1/48) * (a/delta)**4
print(f"\n[1] Two independent loss channels (untwisted bundle, d=0.14mm):")
print(f"  strand-internal skin excess  : {skin_excess:.5f}   -> R/Rdc = {1+skin_excess:.5f}")
print(f"  (matches isolated-strand COMSOL/Bessel value 1.001048)")

# 串扰失衡: 中心-外层互感差 (equal-current circular layout, ln-kernel)
pos = [(0.0, 0.0)]
for k in range(1, 11):
    for j in range(6 * k):
        th = 2 * math.pi * j / (6 * k)
        pos.append((k * p * math.cos(th), k * p * math.sin(th)))
Is = I / N

def sum_ln(px, py, exclude):
    s = 0.0
    for idx, (x, y) in enumerate(pos):
        if idx == exclude:
            continue
        s += math.log(math.hypot(x - px, y - py))
    return s

S_c = sum_ln(0, 0, 0)
eidx = 1 + 6 * 9
ex, ey = pos[eidx]
S_o = sum_ln(ex, ey, eidx)
dL = (mu0 / (2 * math.pi)) * Is * (S_o - S_c)
dX = w * dL
print(f"  center-vs-outer dL = {dL*1e6:.3f} uH/m ; dX = {dX:.3f} ohm/m")
print(f"  dX / R_strand = {dX/R_strand:.3f}  -> 频率越高、d 越细, 该比值越大 (sqrt(f)*d^2 scaling)")
print(f"  observed eta_I = 338% (untwisted circular); this reactance imbalance IS the crosstalk")

# ---------------- [2] 均匀化 Bessel 模型 ----------------
def bessel_R_over_Rdc(a_eff, d_eff):
    """Solid round wire R_ac/R_dc via RK4 integration of Bessel ODE."""
    k = (1 - 1j) / d_eff
    n = 4000
    h = a_eff / n
    r = h / 2
    f = 1 - (k * r) ** 2 / 4
    g = -k * k * r / 2
    for _ in range(n):
        k1f, k1g = g, -g / max(r, 1e-15) - k * k * f
        k2f, k2g = g + h/2*k1g, -(g + h/2*k1g)/max(r + h/2, 1e-15) - k*k*(f + h/2*k1f)
        k3f, k3g = g + h/2*k2g, -(g + h/2*k2g)/max(r + h/2, 1e-15) - k*k*(f + h/2*k2f)
        k4f, k4g = g + h*k3g, -(g + h*k3g)/max(r + h, 1e-15) - k*k*(f + h*k3f)
        f += h/6 * (k1f + 2*k2f + 2*k3f + k4f)
        g += h/6 * (k1g + 2*k2g + 2*k3g + k4g)
        r += h
    J0a, J1a = f, -g / k
    z = (k * a_eff / 2) * J0a / J1a
    return z.real

def J0mag(r, d_eff):
    k = (1 - 1j) / d_eff
    n = 2000
    h = r / n
    rr = h / 2
    f = 1 - (k * rr) ** 2 / 4
    g = -k * k * rr / 2
    for _ in range(n):
        k1f, k1g = g, -g / max(rr, 1e-15) - k * k * f
        k2f, k2g = g + h/2*k1g, -(g + h/2*k1g)/max(rr + h/2, 1e-15) - k*k*(f + h/2*k1f)
        k3f, k3g = g + h/2*k2g, -(g + h/2*k2g)/max(rr + h/2, 1e-15) - k*k*(f + h/2*k2f)
        k4f, k4g = g + h*k3g, -(g + h*k3g)/max(rr + h, 1e-15) - k*k*(f + h*k3f)
        f += h/6 * (k1f + 2*k2f + 2*k3f + k4f)
        g += h/6 * (k1g + 2*k2g + 2*k3g + k4g)
        rr += h
    return abs(f)

d_eff = delta / math.sqrt(fill)
r_hom = bessel_R_over_Rdc(R_env, d_eff)
print(f"\n[2] Homogenized-bundle Bessel model (untwisted):")
print(f"  a_eff/d_eff = {R_env/d_eff:.2f} ; predicted R_ac/R_dc = {r_hom:.3f}")
print(f"  observed: circular 4.403 / hex 4.658 / equal-area solid 4.570")
print(f"  -> 未换位线束 = '准实心导线', 束级趋肤效应主导, 位置电流分串扰是其实股级表现")

rings = [(0, 1)] + [(k, 6 * k) for k in range(1, 11)]
tot = 0.0
wts = []
for kk, n_ in rings:
    r = kk * p
    wt = 1.0 if kk == 0 else J0mag(r, d_eff)
    wts.append((kk, n_, wt))
    tot += n_ * wt
mean = tot / N
obs = [0.0089, 0.0093, 0.0134, 0.0246, 0.0486, 0.0988, 0.2052, 0.4325, 0.9213, 1.9791, 4.2792]
print(f"  ring-by-ring I/mean (predicted | J0 | vs observed COMSOL):")
print(f"  {'k':>3} {'n':>4} {'pred':>8} {'obs':>8}")
for (kk, n_, wt), o in zip(wts, obs):
    print(f"  {kk:3d} {n_:4d} {wt/mean:8.4f} {o:8.4f}")
print(f"  exp decay per ring pitch: {math.exp(-p/d_eff):.3f}")

# ---------------- [3] 完美换位地板 ----------------
# 股内涡流损耗推导(半径 a 圆柱, 横向时谐场 B_rms):
#   E_phi(r) = -(dB/dt) r/2  ->  P_bar = sigma*pi*a^4*w^2*B_rms^2/8  (单位长度)
C = sigma * math.pi * a ** 4 * w ** 2 / 8   # W/m per strand per (T_rms)^2
print(f"\n[3] Perfect-transposition floor (intra-strand proximity only):")
print(f"  per-strand loss coeff C = {C:.0f} W/(m·T_rms^2)")

def disk_floor(R):
    B2 = (mu0 * I) ** 2 / (8 * math.pi ** 2 * R ** 2)
    return N * C * B2 / Pdc

def annulus_floor(R, r1):
    num = den = 0.0
    M = 2000
    for i in range(M):
        r = r1 + (R - r1) * (i + 0.5) / M
        B = mu0 * I * (r * r - r1 * r1) / (2 * math.pi * r * (R * R - r1 * r1))
        num += B * B * r * (R - r1) / M
        den += r * (R - r1) / M
    return N * C * (num / den) / Pdc

fl = disk_floor(R_env)
print(f"  filled disk R=1.575mm: floor = 1 + {fl:.3f} = {1+fl:.3f}")
print(f"    (closed form fill*N*(d/delta)^4/256 = {fill*N*(d/delta)**4/256:.3f}; ring-sum 0.682)")

# validation: untwisted ACTUAL (imbalanced) ring currents -> predicted intra loss vs COMSOL 0.396 W/m
obs_ratio = [0.0089, 0.0093, 0.0134, 0.0246, 0.0486, 0.0988, 0.2052, 0.4325, 0.9213, 1.9791, 4.2792]
B2w = 0.0
for kk in range(11):
    if kk == 0:
        B = 0.0
    else:
        r = kk * p
        Iin = sum(rings[j][1] * obs_ratio[j] for j in range(kk)) * Is
        Iown = rings[kk][1] * obs_ratio[kk] * Is
        B = mu0 / (2 * math.pi * r) * (Iin + Iown / 2)
    B2w += rings[kk][1] * B * B
P_intra_un = C * B2w
print(f"  untwisted intra check: strand-weighted <B_rms^2> = {B2w/N*1e6:.2f} mT^2 -> "
      f"P_intra = {P_intra_un:.3f} W/m (COMSOL observed 0.396)")

print(f"  annulus layouts (A_Cu fixed, fill_ann=0.60):")
for Rmm in [2.0, 2.5, 3.0, 4.0]:
    R = Rmm * 1e-3
    r1 = math.sqrt(max(R * R - Acu / (0.60 * math.pi), 0))
    af = annulus_floor(R, r1)
    print(f"    R={Rmm:4.1f}mm, r1={r1*1e3:5.2f}mm, thickness={(R-r1)*1e3:4.2f}mm "
          f"({(R-r1)/p:.1f} strand pitches): floor = {1+af:.3f}")
print(f"  -> 空心化/扩径可突破实心排布地板, 但 R->inf 时 floor->1 (度量退化, 需直径约束)")

print(f"\n  strand-diameter scan (filled disk, A_Cu fixed, B field RMS):")
print(f"  {'d/mm':>7} {'N':>6} {'floor':>7}")
for dd in [0.05, 0.07, 0.10, 0.14]:
    dmm = dd * 1e-3
    Nd = Acu / (math.pi * (dmm/2) ** 2)
    Cdd = sigma * math.pi * (dmm/2) ** 4 * w ** 2 / 8
    flr = Nd * Cdd * (mu0 * I) ** 2 / (8 * math.pi ** 2 * R_env ** 2) / Pdc
    print(f"  {dd:7.2f} {Nd:6.0f} {1+flr:7.3f}")
print(f"  -> floor ~ 1 + fill*N*(d/delta)^4/256 ~ 1 + (A_Cu/(64*pi*R^2)) * (d/delta)^4 缩细 d 是双刃剑")

# ---------------- [4] 两级束的尺度比较 ----------------
R_sub = R_strand / 27
I_sub = I / 12.15
dX_L2 = w * mu0 * I / (4 * math.pi)
dX_sub = w * mu0 * I_sub / (4 * math.pi)
print(f"\n[4] Two-level bundling scale check:")
print(f"  bundle-level: dX = {dX_L2:.2f} ohm/m vs R_sub(27 strands) = {R_sub:.4f} ohm/m -> ratio {dX_L2/R_sub:.0f}")
print(f"  (若子束内部不换位, 子束-level 串扰同样毁灭性 -> 多级束每一级都必须换位)")
print(f"  within sub-bundle: dX = {dX_sub:.3f} ohm/m vs R_strand = {R_strand:.3f} ohm/m -> ratio {dX_sub/R_strand:.3f}")

# ---------------- [5] 汇总 ----------------
print(f"\n[5] Summary for design-variable selection:")
print(f"  untwisted:    R/Rdc = 4.40 (circular)       [bundle-scale skin + position crosstalk]")
print(f"  perfect transp: floor = 1.68 (filled R=1.575) [intra-strand proximity]")
print(f"  gap 4.40 -> 1.68: 换位可吃掉的收益空间, 全部来自股间串扰均化")
print(f"  beyond 1.68: 需改 d / N / 排布几何 (annulus, d-scan 见上表)")
print(f"  => Q3 设计变量分层: [换位拓扑] 决定 4.40->1.68, [d,N,截面形状] 决定 1.68->更低")