# 現行研究方針：一般型三次元超曲面のLagrangian torus fibration
更新：2026-11-21 / v1.0（旧統合研究から独立・三次元局所稿11、全体稿13）

## 目標・達成基準・優先順位
一般型の滑らかな三次元超曲面 \(X_d\subset\mathbb P^4\), \(d\ge6\) の指定Kähler形式に対し、**全正則fiberがLagrangian \(T^3\) である全域fibration**を構成する。基底はHausdorffかつ境界なし、特異fiberは連結とする。最初の完成目標でspecial Lagrangian性や特異集合の最小余次元は要求しない。三次元の大域存在定理は**未完成**である。

1. **最優先**：構成済みfocus blockを固定して得たnormal ball補完のelliptic face近傍を変更し、境界なしの基底へ接続する。
2. そのouter toric collarから、\(V_0=X_d\setminus\operatorname{int}N(C')\) 側のfibrationへ相対接着する。Whole fibers、作用周期、coorientation、monodromy、Maslov data、curve relationを追跡する。
3. 閉曲線 \(C'\) を有境界曲面 \(C^\circ\) へ切り分けたときの側面 \(\partial C^\circ\times E\) を構成する。相対 \(H_4\) の消滅をWeinstein構造やLagrangian接着の存在定理と誤認しない。
4. 局所模型論文の公表可能性、証明の依存関係、先行研究との差、数値／整数検証と幾何学的証明の区別を精査する。

## 現行資料と正本
- **主たる公表用原稿**：`papers/threefold-lagrangian-torus-fibration.tex`。三次元大域存在問題を扱う本文付き主原稿（全体稿第13稿、以前の組版記録211頁）。局所接着未完了なので完成した大域定理として扱わない。
- **独立した局所模型論文**：`papers/threefold-blocks.tex`（局所稿第11版、以前の組版記録190頁）。Normal-trivial curve、focus block、normal補完、補集合と相対円板群の計算が中心。
- **原稿を支える全分割節・図・検証素材**：`papers/sections/threefold/`、`verification/run_threefold_checks.py`、`verification/verify_focus_exterior_topology.py`、`verification/verify_focus_relative_models.py`、`research/threefold-validation.json`。
- **研究ノート**：`research/notes.md`。現行の定理・仮定・障害・一次文献・再現範囲を整理し、論文の証明を複製しない。
- 複素次元2の構成・五次整数接着・KE適合問題は [surface-lagrangian-torus-fibration](https://github.com/uedakazushi/surface-lagrangian-torus-fibration) を正本とし、必要時は具体的な定理の仮定・出典を引用する。全体の入口は [research-hub](https://github.com/uedakazushi/research-hub)。

## 三次元大域存在定理の先行研究上の位置
最優先の目標は、滑らかな一般型三次元超曲面 \(X_d\subset\mathbf P^4\), \(d\ge6\), 指定Kähler形式に対して、全正則 \(T^3\) fiberがLagrangianであり、連結な特異fiberとHausdorffな境界のない基底を持ち、全体を覆う大域fibrationを構成すること。最初の定理ではspecial Lagrangian、細いdiscriminant、全特異点の非退化性を必須とせず、適切なrelative local modelsを柔軟に選ぶ。

Mikhalkin（Topology 43 (2004), Theorem 1）は**任意のdegree、dimensionのsmooth hypersurfaceにstratified topological torus fibration**を既に構成した。ただしLagrangianとするのは各maximal cellの選ばれた少なくとも一点のfiberであり、全regular fibersではない。Leung–Wan（JDG 83 (2009), Theorem 3.3）はその選ばれたLagrangian fibersの漸近的special性を一般型hypersurfaceで確立。既知のtopological torus fibrationや一般型のLagrangian toriの存在を本研究の新規性として主張しない。

Evans–Mauri（2021）はSNC/negative vertexの局所模型、Abouzaid–Auroux–Katzarkov（2016）はaffine general type hypersurfacesへのmirror/fibration構成、Fernández de Bobadilla–Pełka（2026, CPAM, DOI 10.1002/cpa.70060）は一般Kähler degenerationのA'Campo/部分領域Lagrangian fibrationを与えるが、いずれも所定のcompact general type X_dの全域大域構成を今回確認した範囲で結論しない。Mak–Matessi–Ruddat–Zharkovによるcompact toric CY hypersurfaceの大域構成は2026年11月の公的講演要旨でも発表されているが、対象はCYである。具体的な定理・書誌・適用範囲は `research/notes.md` と三次元全体稿のIntroduction／`sec:literature` に記す。新規性の網羅的確認とCY構成の全証明との比較は引き続き必要。

## 三次元actual focus補集合の境界・第二ホモトピー・巡回被覆

正本は `papers/sections/threefold/local_blocks/focus_exterior_topology.tex`。構成済みの具体的な積埋込 U=C'×K_c⊂N=C'×D⁴⊂X_d を使う。K_cと同じ抽象的homotopy型を持つ任意の埋込へは拡張しない。e=d−2、g=(d−3)(d−4)/2、V=cl(X_d∖U)、V₀=cl(X_d∖N) とする。

- `tnc:focus-exterior-hopf` / `tnc:focus-exterior-homotopy`：反対位相S¹の商 h=(π(|X|²−|Y|²),XY) で、normal exterior E=D⁴∖int K_c を、原点を含むunknotted solid torusの三次元補集合上のHopf circle bundleとして記述する。E≃S¹∨S²∨S³で、outer S³の包含も同定する。
- `tnc:actual-complement-pi1`：V₀の単連結性とvan Kampenからπ₁(V)=Z。境界π₁(C')×Heisenberg群から、曲線loopsとnormal fiberの二つのloops a,bはVで消え、focus基底境界のsection lift tが生成元となる。UとVへの同時写像のkernelはcentral vanishing-cycle subgroup〈a〉。
- `tnc:all-boundary-homology-maps`：両側への全整数homology写像をmarked basisで決定。Regular torus FはU側に残りVで消える。aをt沿いにtransportしたtorus TはUで消えV側に残る。H₂(∂U)→H₂(V)のcokernelはZ/e、normal face [∂K_c]はH₃(V)の位数eのtorsionを生成し、H₃(V)/im boundary≅H₃(X_d)、H₄ boundary→Vは全射。
- `tnc:equivariant-pi2`：Λ=Z[τ,τ⁻¹]をπ₁(V)の群環とすると、π₂(V)≅Λ^(2g+1)⊕Λ/(e(τ−1))。分解は一般に非標準的だがmoduleの同型型を決定した。整数torsion部分はΛ/eΛ。Universal coverのMayer–Vietoris列で、曲線のdegreeから関係e(τ−1)h=0が生じる。
- `tnc:cyclic-cover-h2`：任意のm≥1の連結巡回被覆V_mはH₂(V_m;Z)≅Z^{m(2g+1)+1}⊕(Z/e)^{m−1}。六次ではπ₂(V)=Λ⁷⊕Λ/(4(τ−1))、二重被覆のH₂=Z¹⁵⊕Z/4。

証明を独立に相互点検し、一般の被覆公式は明示的な整数基底変換でも証明した。`python3 verification/verify_focus_exterior_topology.py --json` はd=6,…,20、m=1,…,12の境界行列・巡回被覆関係・群環coinvariantsを整数演算で照合する。有限検算と一般証明は別であり、形式証明ではない。Auroux §5.1、Akbulut arXiv:1103.3137v2のIntroduction/Figure 1、Hatcher Theorems 4.32/4.41・Example 4.45を一次資料で照合した。標準nodal/fishtail模型やHopf理論自体の新規性は主張せず、今回の具体的な組合せの包括的新規性は未確認。

## 三次元：normal ballへの相対補完と円板の接着条件

正本は `papers/sections/threefold/local_blocks/focus_relative_filling.tex`（局所稿第37節、全体稿の同名節）。標準symplectic normal ballで次を証明した。

- `ncc:normal-ball-completion` / `ncc:product-completion`：逆位相circleの商座標A=π(|X|²−|Y|²), w=XYに対し、平面diffeomorphism ψ_Aを内側でw−c、outer sphere近傍でidentityとする。H=|ψ_A(w)|²によりĜ=(A,H)を構成し、K_cの開近傍で元のG_cと一致、Ĝ⁻¹(R)=K_c、outer collarで標準toric、全fiber連結を証明する。Curve quotientとの積はregular Lagrangian T³を持ち、既存のfocus–hyperbolic fibersを保持する。
- `ncc:elliptic-face-limitation`：新しいH=0のfiberはregular curve level上でT²、base germは半空間になる。Compact X_dからHausdorff baseへの大域quotientがこのwhole fiberと近傍の葉層を保つなら、そのgermは外側からの接着だけでは消えない。**Boundary-free目標にはelliptic近傍自体の変更が必要。**
- `tnc:normal-round-handle` / `tnc:relative-pushout`：Eはouter S³から一つのround 1-handleで得られる。逆向きの唯一の2-handleはa^k b^±1を殺し、whiskerでbに正規化できる。V≃V₀∪_i(C'×(S¹∨S²))をV₀の下で同定し、iは[C']=e hを表す。CofiberはS¹∨(∨^{2g+1}S²)∨(∨^{2g+1}S³)∨S⁴。V自体がV₀とのwedgeだという意味ではない。
- `tnc:product-relative-handles` / `tnc:no-weinstein-normal-transition`：C'×Eのouter側からのminimal smooth handle数はindex1,2,3,4で(1,2g+1,2g+1,1)。相対H₄=Z、逆向きでは相対H₅=Zなので、二つの全境界成分を保つWeinstein cobordismはどちら向きにも存在しない。Nonexactな構成や切り分けた境界を否定しない。曲線をpunctureするとこのH₄は消えるが側面∂C°×Eの接着が生じる（`tnc:punctured-relative-handles`）。
- `tnc:normal-relative-disk-module`：適切なreverse 2-handle core class Bについてπ₂(E,Y)=ΛB、∂(fB)=f(1)b+f′(1)a。π₂(E)の像は(τ−1)²ΛBであり、一階差分はa-disk、二階差分はsphere class。任意のb-capをΛ-basisとする主張ではない。
- `tnc:relative-disk-central-extension`：M=π₂(V)に対してπ₂(V,∂U)はMによるπ₁(C')×Z²のcentral extensionで、curve disk liftsは∏[A_i,B_i]=e hを満たす。e h≠0なので非分裂。Degree quotientはEuler数eのcircle bundle P_eのπ₁とZ²の積。Normal pairの全作用は決定したが、全三次元pairでの各curve diskへのt作用とsmooth framingは未決定。
- `ncc:normal-maslov-marking`：normal regular torusのμ(a)=0, μ(b)=2。境界でcleanに交わるnormal C²内のLagrangian diskがma+nbを埋めるならn=0, mA=0。従ってb+kaはこの方法では埋まらない。Elliptic collapseのtransverse symplectic diskやcurve方向にも動くdiskは対象外。

Symingtonのnodal trade、Brownのrelative crossed moduleとHatcherのrelative Hurewicz、Cieliebak–EliashbergのWeinstein index制限を一次資料と照合した。Asimovのround handleの出典は公式書誌を確認し、用いる分解は直接証明した。既知の方法そのものの新規性は主張せず、本研究の具体的相対データの包括的新規性は未確認。全三次元の大域fibrationは未完成である。

## 計算・組版・保存の状態
2026-11-21の分割前の作業記録では、相対円板の243 Laurent計算と1,701 deck作用、15 closed曲線および60 punctured曲線の相対Künneth計算を整数で検証し、三次元局所稿190頁・全体稿211頁をLuaLaTeXで組版し代表頁を表示確認した。**本リポジトリへの移管自体は数学的な再検証を意味しない**。実行した検証と組版をその都度別に記録する。
`latexmk papers/threefold-lagrangian-torus-fibration.tex` と `latexmk papers/threefold-blocks.tex` をリポジトリ直下で用いる。LuaLaTeX・LuaTeX-ja・TikZ・graphicx・日本語フォント等に依存し、PDF・SyncTeXは `build/` に出力してGitで追跡しない。
整数検証は `python3 verification/verify_focus_exterior_topology.py --json`、`python3 verification/verify_focus_relative_models.py --json`、`python3 verification/run_threefold_checks.py`。局所模型の図6点は元の画像とその生成Pythonを保持する。形式証明および外部査読は未実施。原理と今回の具体的な相対構成の包括的新規性は未確認。

## 共通研究原則と保存許可（上記の研究固有方針とは独立）
数学的重要性、新規性、理解の深化、応用可能性、正確性、実現可能性を総合して方針を柔軟に更新する。既知の結果・証明済み・追加仮定付き・未完成・予想・計算観察を明示し、証明の所在と未解決を区別する。文献は一次資料の範囲を精査し、原稿に証明があること、今回の監査、計算確認、形式検証を混同しない。原稿は日本語LuaLaTeXで証明を保持し、現在有効な知識と障害をノートへ整理し、旧版・ログを蓄積しない。
必要なソース、図、検証コード、証明書を軽量かつ再現可能に維持する。研究の区切りには承諾なしで通常コミットし `main` にpush、長い作業では15–30分で保存する。他の作業者の変更を保持し、force push・GitHub Actionsの追加／実行はしない。GitHub正本化後はLibraryへ定常的に二重保存・ZIP保存せず、push不能時の緊急保存または別の指示時のみ使用する。Library既存ファイルを勝手に削除しない。対象研究の現行資料保存後、research-hubの到達点・未解決・次の一手・相互依存・再開先を更新し、対象repoとhubの保存成否を別々に確認する。
