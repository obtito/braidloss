# -*- coding: utf-8 -*-
"""
论文插图重绘脚本（问题一 / 问题二）
数据来源：COMSOL/Q1_COMSOL、Q2_Circular_COMSOL、Q2_RegularTwist_COMSOL 的原始导出 CSV
本脚本只做「读原始数据 -> 绘图」，不重算物理量、不外推、不插补。
输出：论文/figs/*.pdf 与 *.png
"""
import os
import json
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ---------------------------------------------------------------- 全局样式
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
Q1 = os.path.join(ROOT, "COMSOL", "Q1_COMSOL")
Q2C = os.path.join(ROOT, "COMSOL", "Q2_Circular_COMSOL")
Q2T = os.path.join(ROOT, "COMSOL", "Q2_RegularTwist_COMSOL")
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "figs"))
os.makedirs(OUT, exist_ok=True)

C_MAIN = "#1f4e79"   # 主色（深蓝）
C_ALT = "#c00000"    # 对比色（中国习惯：红用于强调/上升）
C_GRN = "#2e7d32"    # 第三方案（绿）
C_GRAY = "#7f7f7f"
C_ORG = "#e07b00"
GRID = dict(color="0.85", lw=0.6, ls="-")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "font.size": 10.5,
    "axes.labelsize": 11,
    "axes.titlesize": 11.5,
    "legend.fontsize": 9.5,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "axes.linewidth": 0.9,
    "figure.dpi": 120,
    "savefig.dpi": 400,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
})


def style(ax, xlabel=None, ylabel=None, title=None):
    ax.grid(True, which="major", **GRID)
    ax.set_axisbelow(True)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="minor", length=2)
    ax.tick_params(which="both", direction="out")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, pad=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("  ->", name)


def load_csv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    return rows


def col(rows, key, dtype=float):
    return np.array([dtype(r[key]) for r in rows])


# ================================================================ 图 1
# 问题一：径向电流密度分布，有限元 vs Bessel 精确解
def fig_q1_radial():
    rows = load_csv(os.path.join(Q1, "radial_theory_comparison.csv"))
    r_mm = col(rows, "r_m") * 1e3
    fem_re = col(rows, "Jz_real_FEM") / 1e6   # A/m^2 -> A/mm^2
    fem_im = col(rows, "Jz_imag_FEM") / 1e6
    bes_re = col(rows, "Jz_real_Bessel") / 1e6
    bes_im = col(rows, "Jz_imag_Bessel") / 1e6
    fem = np.hypot(fem_re, fem_im)
    bes = np.hypot(bes_re, bes_im)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ax = axes[0]
    ax.plot(r_mm, bes, "-", color=C_ALT, lw=3.2, alpha=0.45,
            label="圆柱精确解（Bessel）")
    ax.plot(r_mm, fem, "--", color=C_MAIN, lw=1.5, label="COMSOL 有限元")
    ax.axvline(1.0 - 0.14777165489999476, color=C_GRAY, lw=1.1, ls=":")
    ax.annotate("$r=a-\\delta$\n(147.77 μm)", xy=(1.0 - 0.14777, 6.0),
                xytext=(0.30, 14.0), fontsize=9, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.45", lw=0.9))
    style(ax, xlabel="径向坐标 $r$ / mm",
          ylabel="轴向电流密度有效值 $|J_z|$ / (A·mm$^{-2}$)",
          title="(a) 电流密度的径向分布")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 34)
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    ax.plot(r_mm, fem_re, "-", color=C_MAIN, lw=1.6, label="$J_z$ 实部（FEM）")
    ax.plot(r_mm, fem_im, "-", color=C_ORG, lw=1.6, label="$J_z$ 虚部（FEM）")
    ax.plot(r_mm, bes_re, ":", color=C_MAIN, lw=2.6, alpha=0.5,
            label="$J_z$ 实部（Bessel）")
    ax.plot(r_mm, bes_im, ":", color=C_ORG, lw=2.6, alpha=0.5,
            label="$J_z$ 虚部（Bessel）")
    style(ax, xlabel="径向坐标 $r$ / mm",
          ylabel="电流密度相量分量 / (A·mm$^{-2}$)",
          title="(b) 复数相量分量对照")
    ax.set_xlim(0, 1.0)
    ax.legend(loc="upper left", frameon=False, ncol=1)

    fig.tight_layout()
    save(fig, "q1_radial_J")


# ================================================================ 图 2
# 问题一：网格收敛性
def fig_q1_mesh():
    rows = load_csv(os.path.join(Q1, "mesh_convergence.csv"))
    n = col(rows, "radial_elements", int)
    h = col(rows, "radial_h_um")
    rac = col(rows, "Rac_ohm") * 1e3
    # 解析解参考值（来自 results_summary.json）
    with open(os.path.join(Q1, "results_summary.json"), encoding="utf-8") as f:
        sm = json.load(f)
    rac_ex = sm["Rac_exact_mohm"]
    per_delta = 147.77165489999476 / h

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ax = axes[0]
    ax.axhline(rac_ex, color=C_ALT, lw=1.4, ls="--",
               label=f"解析解 {rac_ex:.6f} mΩ")
    ax.semilogx(n, rac, "o-", color=C_MAIN, lw=1.6, ms=6,
                mfc="white", mew=1.6, label="COMSOL 有限元")
    for xi, yi, hi in zip(n, rac, per_delta):
        ax.annotate(f"{yi:.4f}", xy=(xi, yi), xytext=(0, 7),
                    textcoords="offset points", ha="center", fontsize=8,
                    color="0.25")
    style(ax, xlabel="径向单元数 $N_r$",
          ylabel="交流电阻 $R_{ac}$ / mΩ·m$^{-1}$",
          title="(a) 电阻随网格的收敛")
    ax.set_xlim(6.5, 170)
    ax.set_ylim(19.990, 20.022)
    ax.legend(loc="lower right", frameon=False)

    ax = axes[1]
    err = np.abs(rac - rac_ex) / rac_ex * 100
    ax.loglog(h, err, "s-", color=C_MAIN, lw=1.6, ms=6, mfc="white",
              mew=1.6, label="相对解析解误差")
    ax.loglog(h, err, ":", color=C_GRAY, lw=1.0)
    for xi, yi in zip(h, err):
        ax.annotate(f"{yi:.2e}", xy=(xi, yi), xytext=(4, -11),
                    textcoords="offset points", fontsize=8, color="0.25")
    ax.axvline(7.8125, color=C_ALT, lw=1.1, ls=":")
    ax.annotate("采用网格\n7.8125 μm\n(≈18.9 单元/δ)", xy=(7.8125, 3e-4),
                xytext=(0.0, 0.0), textcoords="axes fraction",
                xycoords="data", fontsize=8.5, color=C_ALT,
                ha="left", va="center")
    style(ax, xlabel="径向单元尺寸 $h$ / μm",
          ylabel="相对解析解误差 / %",
          title="(b) 误差随网格尺寸的变化")
    ax.legend(loc="lower right", frameon=False)

    fig.tight_layout()
    save(fig, "q1_mesh_convergence")


# ================================================================ 图 3
# 问题一：50 kHz - 1 MHz 频率扫描
def fig_q1_freq():
    rows = load_csv(os.path.join(Q1, "frequency_sweep.csv"))
    f = col(rows, "frequency_Hz") / 1e3
    rac = col(rows, "Rac_ohm") * 1e3
    ratio = col(rows, "Rac_over_Rdc")
    delta = col(rows, "delta_m") * 1e6
    loss = col(rows, "loss_W")

    rows2 = load_csv(os.path.join(Q1, "frequency_theory_comparison.csv"))
    f2 = col(rows2, "frequency_Hz") / 1e3
    rac_b = col(rows2, "Rac_Bessel_ohm") * 1e3

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ax = axes[0]
    ax.plot(f2, rac_b, "-", color=C_ALT, lw=3.2, alpha=0.45,
            label="解析解（Bessel）")
    ax.plot(f, rac, "o--", color=C_MAIN, lw=1.5, ms=6, mfc="white",
            mew=1.6, label="COMSOL 有限元")
    ax.axvline(200, color=C_GRAY, lw=1.1, ls=":")
    ax.annotate("统一工况 200 kHz", xy=(200, 22), xytext=(300, 15),
                fontsize=9, color="0.25",
                arrowprops=dict(arrowstyle="->", color="0.45", lw=0.9))
    style(ax, xlabel="频率 $f$ / kHz",
          ylabel="交流电阻 $R_{ac}$ / mΩ·m$^{-1}$",
          title="(a) 交流电阻的频率特性")
    ax.set_xlim(0, 1050)
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    ax.plot(f, ratio, "s-", color=C_MAIN, lw=1.6, ms=6, mfc="white",
            mew=1.6, label="$R_{ac}/R_{dc}$（有限元）")
    ax.axhline(1.0, color=C_GRN, lw=1.2, ls="--", label="理想值 $R_{ac}/R_{dc}=1$")
    ax.axvline(200, color=C_GRAY, lw=1.1, ls=":")
    ax.annotate(f"200 kHz: {ratio[2]:.4f}", xy=(200, ratio[2]),
                xytext=(280, 2.3), fontsize=9, color=C_MAIN,
                arrowprops=dict(arrowstyle="->", color=C_MAIN, lw=0.9))
    style(ax, xlabel="频率 $f$ / kHz",
          ylabel="$R_{ac}/R_{dc}$",
          title="(b) 交直流电阻比的频率特性")
    ax.set_xlim(0, 1050)
    ax.set_ylim(0.8, 8.6)
    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout()
    save(fig, "q1_frequency_sweep")

    # 附：趋肤深度随频率变化（供正文引用，简洁单图）
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.loglog(f, delta, "o-", color=C_MAIN, lw=1.6, ms=6, mfc="white",
              mew=1.6, label="趋肤深度 $\\delta$")
    ax.axhline(147.77165489999476, color=C_ALT, lw=1.2, ls="--",
               label="200 kHz 处 $\\delta$=147.77 μm")
    style(ax, xlabel="频率 $f$ / kHz", ylabel="趋肤深度 $\\delta$ / μm",
          title="铜的趋肤深度随频率的变化")
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    save(fig, "q1_skin_depth")


# ================================================================ 图 4
# 问题二：逐股电流分布（圆形简单成束）
def fig_q2_strand_current():
    rows = load_csv(os.path.join(Q2C, "strand_analysis.csv"))
    ring = col(rows, "radial_ring", int)
    r_mm = col(rows, "radius_mm")
    I_mA = col(rows, "I_rms_magnitude_A") * 1e3
    ideal = 20.0 / 331 * 1e3  # 60.42296 mA

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.0),
                             gridspec_kw={"width_ratios": [1.05, 1.0]})

    # (a) 按半径散点：每根股线的净电流
    ax = axes[0]
    sc = ax.scatter(r_mm, I_mA, c=I_mA, cmap="RdYlBu_r", s=22,
                    edgecolors="0.35", linewidths=0.25, zorder=3,
                    vmin=0, vmax=I_mA.max())
    ax.axhline(ideal, color=C_GRN, lw=1.5, ls="--", zorder=4,
               label=f"理想均分 {ideal:.3f} mA")
    style(ax, xlabel="股线中心半径 $r$ / mm",
          ylabel="单股净电流有效值 $I_i$ / mA",
          title="(a) 331 股电流沿径向的分布")
    ax.set_ylim(-10, 290)
    ax.legend(loc="upper left", frameon=False)
    cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("$I_i$ / mA", fontsize=9.5)
    cb.ax.tick_params(labelsize=8.5)
    cb.outline.set_linewidth(0.6)

    # (b) 分圈平均电流柱状 + 理想线
    ax = axes[1]
    ring_stats = load_csv(os.path.join(Q2C, "ring_statistics.csv"))
    rr = col(ring_stats, "radial_ring", int)
    cnt = col(ring_stats, "count", int)
    meanI = col(ring_stats, "mean_Irms_A") * 1e3
    colors = [C_ALT if i == rr.max() else C_MAIN for i in rr]
    bars = ax.bar(rr, meanI, color=colors, edgecolor="white", lw=0.7,
                  width=0.72, zorder=3)
    ax.axhline(ideal, color=C_GRN, lw=1.5, ls="--", zorder=4,
               label=f"理想均分 {ideal:.3f} mA")
    for b, v, c in zip(bars, meanI, cnt):
        ax.annotate(f"{c}股", xy=(b.get_x() + b.get_width() / 2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=8, color="0.3")
    style(ax, xlabel="同心圆环编号 $k$（半径 $r=0.15k$ mm）",
          ylabel="该圈单股平均电流 $\\bar{I}_k$ / mA",
          title="(b) 各同心圆环的平均股电流")
    ax.set_xticks(rr)
    ax.set_ylim(0, 300)
    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout()
    save(fig, "q2_strand_current")

    # 附：分圈铜损占比（半对数）
    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    frac = col(ring_stats, "loss_fraction_percent")
    ax.semilogy(rr, np.maximum(frac, 1e-6), "o-", color=C_MAIN, lw=1.6,
                ms=6, mfc="white", mew=1.6)
    ax.semilogy(rr[-1], frac[-1], "o", color=C_ALT, ms=9, zorder=5)
    ax.annotate(f"外圈 60 股\n占总铜损 {frac[-1]:.2f}%",
                xy=(rr[-1], frac[-1]), xytext=(-14, -46),
                textcoords="offset points", fontsize=9, color=C_ALT,
                ha="right",
                arrowprops=dict(arrowstyle="->", color=C_ALT, lw=1.0))
    style(ax, xlabel="同心圆环编号 $k$",
          ylabel="该圈铜损占总铜损比例 / %",
          title="铜损沿径向的高度集中")
    ax.set_xticks(rr)
    ax.set_ylim(1e-6, 3e2)
    fig.tight_layout()
    save(fig, "q2_ring_loss")


# ================================================================ 图 5
# 问题二：三方案等铜截面对照 + 频率响应
def fig_q2_compare():
    rows = load_csv(os.path.join(Q2C, "shape_comparison.csv"))
    name_map = {
        "equal_area_solid": "等铜截面实心线",
        "hexagonal_parallel_bundle": "六边形简单成束",
        "circular_parallel_bundle": "圆形简单成束",
    }
    order = ["equal_area_solid", "hexagonal_parallel_bundle",
             "circular_parallel_bundle"]
    by = {r["configuration"]: r for r in rows}
    labels = [name_map[k] for k in order]
    rac = np.array([float(by[k]["Rac_mohm_per_m"]) for k in order])
    loss = np.array([float(by[k]["loss_W_per_m"]) for k in order])
    cols = [C_GRAY, C_ORG, C_MAIN]

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))

    ax = axes[0]
    b = ax.bar(labels, rac, color=cols, edgecolor="white", lw=0.8,
               width=0.6, zorder=3)
    for bi, v in zip(b, rac):
        ax.annotate(f"{v:.4f}", xy=(bi.get_x() + bi.get_width() / 2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=9.5, color="0.2")
    style(ax, ylabel="交流电阻 $R_{ac}$ / mΩ·m$^{-1}$",
          title="(a) 三种结构的交流电阻")
    ax.set_ylim(14.0, 16.2)
    ax.tick_params(axis="x", labelsize=9)

    ax = axes[1]
    # 圆形 vs 六边形 频率响应（真实 COMSOL 数据）
    fc = load_csv(os.path.join(Q2C, "frequency_sweep.csv"))
    fh = load_csv(os.path.join(Q2C, "reference_hexagonal",
                               "frequency_sweep.csv"))
    fc_x = col(fc, "frequency_Hz") / 1e3
    fc_y = col(fc, "Rac_ohm") * 1e3
    fh_y = col(fh, "Rac_ohm") * 1e3
    ax.plot(fc_x, fc_y, "o-", color=C_MAIN, lw=1.6, ms=6, mfc="white",
            mew=1.6, label="圆形简单成束")
    ax.plot(fc_x, fh_y, "s--", color=C_ORG, lw=1.6, ms=6, mfc="white",
            mew=1.6, label="六边形简单成束")
    ax.axvline(200, color=C_GRAY, lw=1.1, ls=":")
    style(ax, xlabel="频率 $f$ / kHz",
          ylabel="交流电阻 $R_{ac}$ / mΩ·m$^{-1}$",
          title="(b) 两种排布方式的频率响应")
    ax.set_xlim(0, 1050)
    ax.legend(loc="upper left", frameon=False)

    fig.tight_layout()
    save(fig, "q2_structure_compare")


# ================================================================ 图 6
# 问题二：规则绞合的损耗分解与三方案总览（含绞合）
def fig_q2_twist():
    # (a) 三方案总览（含规则绞合）
    comp = load_csv(os.path.join(Q2T, "comparison.csv"))
    labels = [r["model"] for r in comp]
    rac = col(comp, "Rac_mohm_per_m")
    loss = col(comp, "copper_loss_W_per_m")
    ratio = col(comp, "Rac_over_Rdc")
    cols = [C_GRAY, C_ORG, C_MAIN]

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))

    ax = axes[0]
    b = ax.bar(labels, rac, color=cols, edgecolor="white", lw=0.8,
               width=0.6, zorder=3)
    for bi, v in zip(b, rac):
        ax.annotate(f"{v:.4f}", xy=(bi.get_x() + bi.get_width() / 2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=9.5, color="0.2")
    style(ax, ylabel="交流电阻 $R_{ac}$ / mΩ·m$^{-1}$",
          title="(a) 三方案交流电阻对照")
    ax.set_ylim(13.9, 15.8)
    ax.tick_params(axis="x", labelsize=8.6)

    # (b) 圆形成束铜损三分量分解
    ax = axes[1]
    with open(os.path.join(Q2C, "results_summary.json"), encoding="utf-8") as f:
        sc = json.load(f)
    parts = [sc["ideal_dc_loss_W_per_m"],
             sc["strand_imbalance_excess_loss_W_per_m"],
             sc["within_strand_excess_loss_W_per_m"]]
    pnames = ["理想直流基准\n$I^2R_{dc}$",
              "股间电流不均\n（幅值+相位）",
              "股内电流分布\n不均"]
    pcols = [C_GRN, C_ALT, C_ORG]
    bottom = 0.0
    for v, nm, cc in zip(parts, pnames, pcols):
        ax.bar(["圆形简单成束"], [v], bottom=bottom, color=cc,
               edgecolor="white", lw=0.8, width=0.42, label=nm, zorder=3)
        ax.annotate(f"{v:.4f} W/m\n({v/sum(parts)*100:.2f}%)",
                    xy=(0, bottom + v / 2), ha="center", va="center",
                    fontsize=8.6, color="white", fontweight="bold")
        bottom += v
    ax.annotate(f"合计 {sum(parts):.4f} W/m", xy=(0, bottom),
                xytext=(0, 6), textcoords="offset points",
                ha="center", fontsize=9.5, color="0.2")
    style(ax, ylabel="每米铜损 / W·m$^{-1}$",
          title="(b) 圆形成束铜损的构成分解")
    ax.set_ylim(0, 7.4)
    ax.legend(loc="upper right", frameon=False, fontsize=8.4)

    fig.tight_layout()
    save(fig, "q2_twist_and_decomposition")


# ================================================================ 图 7
# 问题二：规则绞合的分层电流与损耗
def fig_q2_twist_ring():
    rows = load_csv(os.path.join(Q2T, "ring_statistics.csv"))
    rr = col(rows, "ring", int)
    r_mm = col(rows, "radius_mm")
    cnt = col(rows, "count", int)
    meanI = col(rows, "current_mean_rms_mA")
    share = col(rows, "copper_loss_share_percent")

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 3.9))

    ax = axes[0]
    ax.bar(rr, meanI, color=C_MAIN, edgecolor="white", lw=0.7,
           width=0.72, zorder=3)
    ax.axhline(20.0 / 331 * 1e3, color=C_GRN, lw=1.5, ls="--", zorder=4,
               label="理想均分 60.423 mA")
    for b, v, c in zip(ax.patches, meanI, cnt):
        ax.annotate(f"{c}股", xy=(b.get_x() + b.get_width() / 2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", fontsize=8, color="0.3")
    style(ax, xlabel="径向层编号 $k$（半径 $r=0.15k$ mm）",
          ylabel="该层单股平均电流 $\\bar{I}_k$ / mA",
          title="(a) 规则绞合：分层平均电流")
    ax.set_xticks(rr)
    ax.set_ylim(0, 290)
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    ax.semilogy(rr, np.maximum(share, 1e-4), "o-", color=C_MAIN, lw=1.6,
                ms=6, mfc="white", mew=1.6)
    ax.semilogy(rr[-1], share[-1], "o", color=C_ALT, ms=9, zorder=5)
    ax.annotate(f"外圈 60 股\n占总铜损 {share[-1]:.2f}%",
                xy=(rr[-1], share[-1]), xytext=(-14, -46),
                textcoords="offset points", fontsize=9, color=C_ALT,
                ha="right",
                arrowprops=dict(arrowstyle="->", color=C_ALT, lw=1.0))
    style(ax, xlabel="径向层编号 $k$",
          ylabel="该层铜损占总铜损比例 / %",
          title="(b) 规则绞合：铜损的径向集中")
    ax.set_xticks(rr)
    ax.set_ylim(1e-4, 3e2)

    fig.tight_layout()
    save(fig, "q2_twist_ring")


# ================================================================ 图 8
# 问题二：绞合 vs 成束 的逐股电流差异（按半径对照）
def fig_q2_strand_compare():
    rc = load_csv(os.path.join(Q2C, "strand_analysis.csv"))
    rt = load_csv(os.path.join(Q2T, "strand_currents_200kHz.csv"))
    r_c = col(rc, "radius_mm")
    I_c = col(rc, "I_rms_magnitude_A") * 1e3
    r_t = col(rt, "radius_mm")
    I_t = np.hypot(col(rt, "I_real_rms_A"), col(rt, "I_imag_rms_A")) * 1e3

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.scatter(r_c, I_c, s=20, color=C_ORG, alpha=0.75,
               edgecolors="none", label="圆形简单成束（直线）")
    ax.scatter(r_t, I_t, s=20, color=C_MAIN, alpha=0.75, marker="^",
               edgecolors="none", label="规则绞合（绞距 40 mm）")
    ax.axhline(20.0 / 331 * 1e3, color=C_GRN, lw=1.5, ls="--",
               label="理想均分 60.423 mA")
    style(ax, xlabel="股线中心半径 $r$ / mm",
          ylabel="单股净电流有效值 $I_i$ / mA",
          title="两种结构下股电流沿径向的分布对照")
    ax.set_ylim(-10, 290)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    save(fig, "q2_strand_compare")


if __name__ == "__main__":
    print("输出目录：", OUT)
    print("问题一插图：")
    fig_q1_radial()
    fig_q1_mesh()
    fig_q1_freq()
    print("问题二插图：")
    fig_q2_strand_current()
    fig_q2_compare()
    fig_q2_twist()
    fig_q2_twist_ring()
    fig_q2_strand_compare()
    print("完成。")
