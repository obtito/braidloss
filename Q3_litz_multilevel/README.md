# 问题三成果：传统多级递归绞合

2026-09-27 续作完成。入口：[简报](进度与控制变量.md) · [完整 Markdown 报告](问题三_结果与控制变量.md) · [PDF 验证报告](问题三_续作验证报告.pdf) · [最终数据](final_results.json)。

固定 64 股、6 mm²、200 kHz、20 A RMS。推荐参数为 4×16 两级绞合，内级 -16 mm、外级 +128 mm。改进幅度与精度边界以报告和最终数据为准，未证明唯一或全局最优。

![实际场对比](figures/final_comparison.png)

- 正式冻结基准：`data/raw/baseline_64_av_m3/`。
- 本次最终推荐：`data/raw/recommended_m325/`。
- [129 组筛选表](parameter_search.csv)、[方法说明](MODEL_METHOD.md)、[预定范围](DESIGN_SPACE.md)。
- [后续验证限制](data/rom_post_confirmation.json)：新三级样本超过原定 5% ROM 门槛，保留历史筛选、限制未经新验证的继续预测。
- 数据复核结果保存为 `delivery_audit.json`；核查入口 `src/audit_delivery.py`。

复现：从本目录运行 `python src/finalize_study.py`、`python src/write_report.py`、`python src/build_report_pdf.py`、`python src/audit_delivery.py`。PDF 脚本需 ReportLab；也可使用本机 Codex 随附的文档 Python。需要新的 COMSOL 计算时采用新案例名，保留既有原始证据。

本目录早期 343 股几何示意、失败试算及外部旧目录 `litz_q3` 的方格换位结论不作为本轮推荐依据。本次成果是问题三既定范围的研究闭环，非整套四题与网站的全部交付。
