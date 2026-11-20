# Lagrangian torus fibrations on general-type threefolds
2026-11-21 / v1.0 — [現行研究方針](00_research_plan.md) · [研究ノート](research/notes.md) · [研究全体の入口](https://github.com/uedakazushi/research-hub)

滑らかな一般型三次元超曲面 \(X_d\subset\mathbb P^4\), \(d\ge6\) の指定Kähler形式に対して、**全正則fiberがLagrangian \(T^3\) の大域fibration**を構成する研究。**大域存在問題は未完成**。既知のMikhalkinのstratified topological torus fibrationと、今回目標とする全正則fiberのLagrangian性を区別する。

## 公表用原稿
- [三次元の主原稿](papers/threefold-lagrangian-torus-fibration.tex)（全体稿第13稿）。Boundary-freeな大域構成を目標とし、現時点では局所模型・補集合・相対接着に関する証明と未解決点を含む。
- [局所模型の独立原稿](papers/threefold-blocks.tex)（第11版）。閉normal-trivial symplectic曲線、focus–hyperbolic block、normal ball補完、補集合の整数位相、相対handles・disk groups、適用制約。全ての分割節・原画像・生成コードは [papers/sections/threefold/](papers/sections/threefold/) に保持。
- 原稿は独立に読めるよう整える。新規成果の包括的新規性は未確認。外部査読・形式証明は未実施。

## 現行の到達点と未解決の核心
Curve degree \(e=d-2\)、genus \(g=(d-3)(d-4)/2\) のnormal-trivialな実際のsymplectic曲線と、whole focus–hyperbolic fibers \((d-2)(d-5)\) 個の局所Lagrangian \(T^3\) blockを構成。Normal ball上で既存blockの近傍を固定してouter collarをtoricにする補完を得たが、新しい**elliptic base boundary**が残る。相対 \(H_4\)／逆向き \(H_5\) のため、閉曲線全体の標準Weinstein移行はできない。Normal相対円板群は \(\pi_2(E,Y)=\mathbb Z[\tau^{\pm1}]B\)、\(\partial(fB)=f(1)b+f'(1)a\)。補集合のequivariant \(\pi_2\)・巡回被覆homologyも決定しているが、**elliptic近傍の変更、管状補集合への接続、大域接着**は未解決。証明の範囲と仮定は [研究方針](00_research_plan.md) と各TeX原稿を参照。

## 再現・ビルド
- `latexmk papers/threefold-lagrangian-torus-fibration.tex`、`latexmk papers/threefold-blocks.tex`（LuaLaTeX、LuaTeX-ja、日本語フォント、graphicx/TikZなど）。生成PDFとSyncTeXは `build/` に出力し、Gitに追跡しない。
- `python3 verification/verify_focus_exterior_topology.py --json`、`python3 verification/verify_focus_relative_models.py --json`、`python3 verification/run_threefold_checks.py`。検証コードと原図、旧作業で得た記録を保持。有限計算は一般のsymplectic存在証明を代替しない。
- 分離前2026-11-21の組版記録：局所稿190頁、全体稿211頁。これは**新リポジトリでの再組版記録ではない**。

## 関連研究
複素次元2の大域Lagrangian \(T^2\) 構成、五次整数アフィン接着、曲面のKEとの適合は [surface-lagrangian-torus-fibration](https://github.com/uedakazushi/surface-lagrangian-torus-fibration) の `main` を正本とする。必要な結果を引用する際は仮定と具体的な定理を明記する。
