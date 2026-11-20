# 三次元の公表用原稿
更新：2026-11-21

- `threefold-lagrangian-torus-fibration.tex`：一般型三次元超曲面全体へのLagrangian torus fibrationを目標とする主原稿（大域存在は未完成）。分割元では `threefold-fibration.tex` と呼んでいた同じ全体稿を、主原稿の命名規則で配置。
- `threefold-blocks.tex`：独立した局所模型論文。Normal-trivial curves、focus blocks、normal補完・相対接着の計算と障害。局所稿には一般型三次元の全域存在定理は含まない。
- `sections/threefold/`：両論文から必要な節、生成図、検証入力。TeX root指示は主原稿を指す。局所稿も独立にlatexmkで組版可能。
- `../00_research_plan.md` および `../research/notes.md` が現行状態と未解決点の入口。

組版はリポジトリルートで `latexmk papers/threefold-lagrangian-torus-fibration.tex` または `latexmk papers/threefold-blocks.tex`。結果は `build/` に出力する。
