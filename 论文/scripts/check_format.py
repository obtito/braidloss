# -*- coding: utf-8 -*-
"""
论文格式自检脚本
================
用途：编译完成后检查 main.tex / main.pdf 的常见格式陷阱。
运行：python scripts/check_format.py
      （需在 论文/ 目录下运行，且 main.log 与 main.pdf 已存在）

检查项
------
A. 源码级（main.tex）
   A1 可选参数 [...] 内出现全角标点（LaTeX 语法区必须半角）
   A2 \\mathrm/\\mathbf/... 内出现非 ASCII（Times New Roman 数学直立体缺字形）
   A3 行内数学模式内出现裸中文（应包 \\text{}）
   A4 宏重复定义
   A5 \\label 重复
   A6 \\ref/\\eqref 指向不存在的 label
   A7 行内 $ 数量为奇数（疑似未闭合）
   A8 \\U{} 内出现希腊字母/非 ASCII（应改用 \\Uohm/\\Umohm/\\Uohmm 等专用宏）
   A9 行尾多余空白
B. 编译日志（main.log）
   B1 错误（以 ! 开头）
   B2 Overfull box（真正需要修排版）
   B3 Missing character（真正的缺字，需换字体或改写法）
C. PDF 输出（main.pdf）
   C1 未解析引用 ??
   C2 关键字符码位核对：°(U+00B0) / Ω(U+03A9)
      —— 若度符号提取为 U+25E6 说明用了 ^\\circ 而非 \\textdegree
      —— 若 Ω 缺失或为空白说明 \\Omega 被误放进 \\mathrm
   C3 私有区字符（PUA U+E000~U+F8FF）：仅出现在矩阵页属正常，
      是 CMEX10 括号件的文本提取伪影，不是渲染缺陷。
"""
import collections
import os
import re
import sys

FW = "，。；：？！（）【】《》“”‘’、％＃＆"
TEX = "main.tex"
LOG = "main.log"
PDF = "main.pdf"


def hr(title):
    print(f"\n{'=' * 62}\n{title}\n{'=' * 62}")


# ---------------------------------------------------------------- A 源码级
def check_source():
    hr("A. 源码级检查 (main.tex)")
    if not os.path.exists(TEX):
        print(f"  !! 找不到 {TEX}")
        return
    lines = open(TEX, encoding="utf-8").read().split("\n")
    issues = []

    for i, ln in enumerate(lines, 1):
        # A1 可选参数内全角标点
        for m in re.finditer(r"\[([^\[\]]*)\]", ln):
            seg = m.group(1)
            if any(ch in FW for ch in seg) and not re.match(r"^\s*[0-9.,\s]*$", seg):
                issues.append((i, "A1 可选参数内全角标点", seg[:60]))
        # A2 \mathrm{} 等内非 ASCII
        for m in re.finditer(r"\\(mathrm|mathbf|mathit|mathsf|mathtt)\{([^}]*)\}", ln):
            bad = [c for c in m.group(2) if ord(c) > 127]
            if bad:
                issues.append((i, f"A2 \\{m.group(1)}{{}} 内非 ASCII", f"{bad}"))
        # A3 数学模式内裸中文
        for m in re.finditer(r"\$([^$]*)\$", ln):
            body = re.sub(r"\\text\{[^}]*\}", "", m.group(1))
            if any("\u4e00" <= c <= "\u9fff" for c in body):
                issues.append((i, "A3 数学模式内裸中文", m.group(1)[:50]))
        # A7 行内 $ 奇数
        s = re.sub(r"\\\$", "", ln)
        if s.count("$") % 2 == 1 and not s.lstrip().startswith("%"):
            issues.append((i, "A7 行内 $ 数量为奇数", ln.strip()[:60]))
        # A8 \U{} 内希腊/非 ASCII
        for m in re.finditer(r"\\U\{([^}]*)\}", ln):
            if re.search(r"[^\x00-\x7F]|Omega|circ", m.group(1)):
                issues.append((i, "A8 \\U{} 内含希腊/非 ASCII", m.group(1)))
        # A9 行尾空白
        if ln != ln.rstrip() and ln.strip():
            issues.append((i, "A9 行尾多余空白", repr(ln[-8:])))

    # A4 宏重复定义
    defs = {}
    for i, ln in enumerate(lines, 1):
        for m in re.finditer(r"\\(?:re)?newcommand\{\\([A-Za-z]+)\}", ln):
            n = m.group(1)
            if n in defs and not ln.lstrip().startswith("%"):
                issues.append((i, "A4 宏重复定义", f"\\{n} 首见于 L{defs[n]}"))
            defs.setdefault(n, i)

    # A5/A6 label
    labels = {}
    for i, ln in enumerate(lines, 1):
        for m in re.finditer(r"\\label\{([^}]*)\}", ln):
            if m.group(1) in labels:
                issues.append((i, "A5 label 重复", m.group(1)))
            labels.setdefault(m.group(1), i)
    for i, ln in enumerate(lines, 1):
        for m in re.finditer(r"\\(?:ref|eqref|autoref)\{([^}]*)\}", ln):
            if m.group(1) not in labels:
                issues.append((i, "A6 引用不存在的 label", m.group(1)))

    if not issues:
        print("  [OK] 源码级未发现问题")
        return
    by = {}
    for ln, kind, d in issues:
        by.setdefault(kind, []).append((ln, d))
    for kind, items in sorted(by.items()):
        print(f"  [{kind}] {len(items)} 处")
        for ln, d in items[:8]:
            print(f"      L{ln}: {d}")
        if len(items) > 8:
            print(f"      ... 另有 {len(items) - 8} 处")


# ---------------------------------------------------------------- B 日志
def check_log():
    hr("B. 编译日志检查 (main.log)")
    if not os.path.exists(LOG):
        print(f"  !! 找不到 {LOG}")
        return
    log = open(LOG, encoding="utf-8", errors="replace").read()

    errs = [l for l in log.split("\n") if l.startswith("!")]
    print(f"  B1 错误            : {len(errs)}")
    for e in errs[:5]:
        print(f"      {e}")

    ovf = re.findall(r"Overfull \\hbox \(([\d.]+)pt too wide\)", log)
    ovf_v = re.findall(r"Overfull \\vbox", log)
    print(f"  B2 Overfull hbox   : {len(ovf)}" + (f"  最宽 {max(map(float, ovf)):.1f}pt" if ovf else ""))
    print(f"     Overfull vbox   : {len(ovf_v)}")

    miss = re.findall(r"Missing character: There is no (.{1,30})", log)
    print(f"  B3 Missing char    : {len(miss)}")
    for m in collections.Counter(miss).most_common(5):
        print(f"      {m[1]}x  {m[0]!r}")

    und = log.count("undefined")
    print(f"  B4 未定义引用      : {und}")


# ---------------------------------------------------------------- C PDF
def check_pdf():
    hr("C. PDF 输出检查 (main.pdf)")
    if not os.path.exists(PDF):
        print(f"  !! 找不到 {PDF}")
        return
    try:
        from pypdf import PdfReader
    except ImportError:
        print("  !! 未安装 pypdf，跳过（pip install pypdf）")
        return

    r = PdfReader(PDF)
    print(f"  页数: {len(r.pages)}")

    qq = [i for i, p in enumerate(r.pages, 1) if "??" in (p.extract_text() or "")]
    print(f"  C1 含未解析引用 ?? 的页: {qq if qq else '无'}")

    cnt = collections.Counter()
    pua = collections.defaultdict(list)
    for i, p in enumerate(r.pages, 1):
        for ch in (p.extract_text() or ""):
            cp = ord(ch)
            if cp in (0x00B0, 0x25E6, 0x2218, 0x03A9):
                cnt[cp] += 1
            if 0xE000 <= cp <= 0xF8FF:
                pua[cp].append(i)

    print("  C2 关键字符:")
    print(f"      ° U+00B0 真度符号 : {cnt[0x00B0]}   <- 期望 >0")
    print(f"      ◦ U+25E6 白子弹   : {cnt[0x25E6]}   <- 若在温度处出现说明误用 ^\\circ")
    print(f"      ∘ U+2218 环算子   : {cnt[0x2218]}")
    print(f"      Ω U+03A9 欧米伽   : {cnt[0x03A9]}   <- 期望 >0（若为 0 说明 \\Omega 被放进 \\mathrm）")

    if pua:
        pages = sorted({p for v in pua.values() for p in v})
        print(f"  C3 私有区字符: {len(pua)} 个码位，出现在页 {pages}")
        print("      （若仅出现在含矩阵的页，属 CMEX10 括号件提取伪影，非渲染缺陷）")
    else:
        print("  C3 私有区字符: 无")


if __name__ == "__main__":
    print("论文格式自检 —— 当前目录:", os.getcwd())
    check_source()
    check_log()
    check_pdf()
    hr("自检完成")
