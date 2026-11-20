# 現行研究ノート：一般型三次元超曲面のLagrangian torus fibration
更新：2026-11-21 / v1.0

大域存在は未解決である。証明の正本は `papers/threefold-lagrangian-torus-fibration.tex` および `papers/threefold-blocks.tex` とその分割節。曲面のKE退化と整数接着は別の [曲面研究](https://github.com/uedakazushi/surface-lagrangian-torus-fibration) を参照する。以下は現行の文献比較、相対構成、重要な障害、整数検証の要点であり、論文の全文証明は重複させない。

## 三次元一般型のLagrangian torus fibration：一次文献比較（2026-11-20）
最重要の先行結果は以下。**位相的stratified torus fibration**、**選んだLagrangian fiber**、**全regular fibersがLagrangianな大域fibration**、**SNC degenerationの部分領域のLagrangian fibration**を混同しない。

- G. Mikhalkin, *Decomposition into pairs-of-pants for complex algebraic hypersurfaces*, Topology 43 (2004), 1035–1065, DOI 10.1016/j.top.2003.11.006, arXiv:math/0205011。**Theorem 1**は全次数・全次元のsmooth hypersurface X_dにstratified topological T^n-fibrationを作り、底はwedge of p_g(X_d) spheres S^nのhomotopy型を持つ。**Lagrangian性は各maximal n-cell内の少なくとも一点のfiberに限る**。特に一般型threefold d≥6のtopological torus decompositionは既知。全regular fibersが所定ωでLagrangianという結論はない。
- N. C. Leung and T. Y. H. Wan, *Calabi-Yau components in general type hypersurfaces*, J. Differential Geom. 83 (2009), 43–74, DOI 10.4310/jdg/1253804351, arXiv:0807.1784。**Theorem 2.4**がMikhalkinの定理を再掲し、**Theorem 3.3**は選ばれたLagrangian fibersが正則n形式に関してasymptotically specialになることを示す。全正則fiberのLagrangian性は証明しない。
- J. D. Evans and M. Mauri, *Constructing local models for Lagrangian torus fibrations*, Ann. Henri Lebesgue 4 (2021), 537–570, DOI 10.5802/ahl.80。Affine SNC complementのlocal model、3-fold negative vertex、one-dimensional stratum近傍。Projective general type X_dのglobal fibrationではない。
- M. Abouzaid, D. Auroux and L. Katzarkov, *Lagrangian fibrations on blowups of toric varieties and mirror symmetry for hypersurfaces*, Publ. Math. IHÉS 123 (2016), 199–282, DOI 10.1007/s10240-016-0081-9。Affine hypersurfaces of general typeへ応用するが、blowup of V×C along H×0のfibrationを、compact H⊂P^4そのものへのglobal fibrationと同一視できない。
- J. Fernández de Bobadilla and T. Pełka, *Fibrations by Lagrangian Tori for Maximal Calabi–Yau Degenerations and Beyond*, Commun. Pure Appl. Math., first published 8 July 2026, article e70060, DOI 10.1002/cpa.70060, arXiv:2312.13248v4。A'Campo boundaryから近傍fiberへsymplectic transportをする一般Kähler degeneration方法。Maximal CYではgeneric regionのLagrangian torus fibrationであり、一般型compact threefoldのglobal coverageには直接至らない。本研究のthreefold-fibration.texにTheorem 6.14の部分適用が既にある。
- C. Y. Mak, D. Matessi, H. Ruddat and I. ZharkovのCalabi–Yau hypersurfaceへのglobal Lagrangian fibrationは2026年1月/6月/11月の公的**講演要旨**で発表。2026年11月12日のMak講演要旨（IMJ-PRG）はiron coefficientsとconvex non-strict potential、Liouville flowを明示している。確認した範囲でcomplete published proof/preprintは未取得。対象はCY/anticanonicalであり一般型には自動適用できない。公式要旨 https://www.imj-prg.fr/gestion/evenement/affEvenement/45 と https://www.mittag-leffler.se/seminar/helge-ruddat-tba/ 。

**目標の新規性の正確な位置**：smooth general type X_d⊂P^4, d≥6, 指定ω=aHについて、全regular T^3 fibersがLagrangianで、明示されたsingular local modelsをもち、全X_dを覆い、境界なし、連結fiber、Hausdorff baseとなるglobal mapを構成する。Mikhalkinのtopological fibrationとLeung–Wanの部分的special性より強い。しかし「初めての一般型torus fibration」「一般型のLagrangian torusが未知られ」などの主張は誤り。新規性の網羅的な保証は未完了。より弱いgeneral singular setを許しても大域存在を優先する。現在未完成の相対filling、curvature-flat normal lifts、reduced surfaceの両端葉層一致、actual coverageの証明を必要とする。

## 三次元：局所稿の参照先
`papers/threefold-blocks.tex`第25・33節：positive normal shell、membrane/focus collars。第26節：smooth focus Cartanとnormal degreeの障害。第34–35節：compact focus blockとnormal-trivial closed curve。第36節：actual補集合の境界写像・equivariant π₂・巡回被覆。第37節：normal ballへの相対補完・relative handlesとdisk group・接着の障害。第38節：weighted collars。Normal exterior Eのfibrationは構成済み。旧48 ports、V₀への延長、elliptic faceを消すboundary-free接着は未完成。全体稿は `papers/threefold-lagrangian-torus-fibration.tex`。

## 三次元：管状補集合とfocus境界の整数位相
**曲線の構成：** d≥6、e=d−2。平面曲線C_eを含むsmooth hypersurfaceはF=q_eR_2+z_3G_{d−1}+z_4H_{d−1}のfirst-jet生成から取れる。Normal完全列より det N_{C_e/Y_d}=O_{C_e}(e+2−d)=O_{C_e}。Ehresmann–Moser移送したC'⊂X_dのsymplectic normalは自明、genusはg=(d−3)(d−4)/2。一般のnormal-trivial symplectic curve（degree e, genus g）には2g−2=(d−5)eが必要で、d偶数ならe偶数（主稿の tnc:degree-genus-parity）。Holomorphic normal rank-two bundleそのものの自明性は主張しない。

**管状補集合の整数homology：** 小さいnormal disk bundle N(C')の補集合V_0=X_d＼int N(C')では、Lefschetz・Thom長完全列により
\[
H_3(V_0;\mathbb Z)\simeq \mathbb Z^{b_d}\oplus\mathbb Z/(d-2),\qquad
H_4(V_0;\mathbb Z)\simeq\mathbb Z^{2g},\qquad
b_d=d^4-5d^3+10d^2-10d+4.
\]
Torsionはnormal meridian S³のclassから生じる（tnc:tubular-homology）。六次ではH_3(V_0)=Z^{520}⊕Z/4、H_4(V_0)=Z^6。X_dの整数cohomologyの一次資料はO. Randal-Williams, Math. Ann. 391 (2025), 1965–2003, §2.1, DOI 10.1007/s00208-024-02951-4。

**Actual focus blockの整数homology（証明済み）：** $K_c$の座標を$w=XY$、$A=P_X-P_Y$とすると、$w$を保持して$A\to0$へ変形収縮できる。$A=0$は$w$-annulusと円周の積で、$w=0$の一円周を潰したものだから
\[
K_c\simeq S^1\vee S^2
\]
（tnc:focus-homotopy）。従って実際のblockの補集合$\mathcal V=\overline{X_d\setminus\operatorname{int}(C'\times K_c)}$に相対Poincaré–Lefschetz双対性を適用すると
\[
H_1(\mathcal V)=\mathbb Z,\quad H_2(\mathcal V)=\mathbb Z^{2g+2},\quad
H_3(\mathcal V)=\mathbb Z^{b_d+2g+1}\oplus\mathbb Z/(d-2),\quad
H_4(\mathcal V)=\mathbb Z^{2g+1}.
\]
六次では$(b_1,b_2,b_3^{\rm free},b_4)=(1,8,527,7)$、$H_3$のtorsionは$\mathbb Z/4$（tnc:actual-complement-homology）。Torsionを表す3-cycleは固定した曲線上の一点での$\partial K_c$。管状補集合$V_0$と$\mathcal V$は同じEuler数だが自由rankが異なる（tnc:complement-distinction）。Ambientの$\pi_3(X_d)$に現れる既知の$\mathbb Z/d$は別の不変量。

**Actual normal face：** $\partial K_c$はfocus点一つのmonodromy $\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)$を持つtorus bundleで、$\pi_1$は整数Heisenberg群、$H_1=H_2=\mathbb Z^2$（tnc:focus-boundary-homology）。$\partial\mathcal U=C'\times\partial K_c$のBetti数は$(1,2g+2,4g+3,4g+3,2g+2,1)$、六次では$(1,8,15,15,8,1)$（tnc:focus-boundary-betti）。

**Euler予算：** χ(U)=−(d−2)(d−5)、χ(V)=−(d⁴−5d³+9d²−3d−10)。六次ではχ(V)=−512。Euler数非零fiberが有限個の孤立基底点にだけ存在してすべて−1なら512個必要（tnc:complement-singularity-budget）。正次元特異strataではcompact-support Euler数で重み付けし、一般のglobal Lagrangian fibrationの実現数には外挿しない。

**境界包含・被覆空間の現在の到達点：** 証明の正本は `papers/sections/threefold/local_blocks/focus_exterior_topology.tex`。Hopf商によってnormal exterior E≃S¹∨S²∨S³を同定し、V=V₀∪_{C'×S³}(C'×E)を用いてπ₁(V)=Zを証明した。境界のcurve loopsとa,bはVで消え、focus基底boundary lift tは生成元。両側への同時π₁写像のkernelはHeisenberg中心〈a〉である。全整数境界表も決定した（`tnc:all-boundary-homology-maps`）。F=regular torusはU側、T=aをt沿いにtransportしたtorusはV側に残る。H₂ boundary→VのcokernelはZ/e、[∂K_c]はH₃(V)の位数eのtorsion、H₃(V)/im boundary=H₃(X_d)、H₄ boundary→Vは全射。Framing・orientationを固定するが、H₃(X_d)のliftは非標準的なのでcanonicalなsplittingと混同しない。

Λ=Z[τ,τ⁻¹]に対し `tnc:equivariant-pi2` はπ₂(V)=Λ^(2g+1)⊕Λ/(e(τ−1))を与える。Universal coverのMVで[C']=e hとなり、deck不変性からe(τ−1)h=0が生じる。残りのaugmentation ideal (τ−1)Λはfreeなのでextensionがsplitする。整数torsion部分はΛ/eΛだが、coinvariants H₂(V)=Z^(2g+2)はfreeである。`tnc:cyclic-cover-h2` は全m≥1についてH₂(V_m)=Z^{m(2g+1)+1}⊕(Z/e)^{m−1}を証明する。有限coverではH₁(C')⊗H₁(E_m)=Z^(2g)の項を落としてはいけない。六次二重被覆はZ¹⁵⊕Z/4。

**再現と先行研究：** `python3 verification/verify_focus_exterior_topology.py --json` は17個のSmith algorithm自己照合、d6–20の60境界行列、m1–12の180被覆関係行列と180群環coinvariant行列を整数演算で確認する。一般の式は原稿の証明による。Auroux arXiv:0706.3207 §5.1のnodal model、Akbulut arXiv:1103.3137v2 Introduction/Figure 1のS⁴のfishtail分解、Hatcher Theorems 4.32/4.41・Example 4.45を一次資料で確認した。既知の模型・道具と今回の具体計算を区別し、包括的新規性は未確認。

## 三次元：normal補完と相対円板の現在の知識

証明は `papers/sections/threefold/local_blocks/focus_relative_filling.tex` に集約する。

**Positive construction。** `ncc:normal-ball-completion` は、K_c⊂int D⁴の余裕 a_*²+4π²(c+√q₊)²<L²を使い、商平面上のcompactly supported flowでψ_Aを作る。Inner blockを含むfilled diskではw−c、outer collarではidentityであり、H=|ψ_A(XY)|²はAとPoisson可換になる。ψ_A自体にsymplectic性を要求しない点が重要。全fiberの連結性とĜ⁻¹(R)=K_cは各reduced sliceのdiffeomorphism性から従う。C′との積は既存focus–hyperbolic fibersを保存する。商座標の図はintegral-affine図ではない。

**大域化の障害。** H=0は真のelliptic faceであり、そのwhole T² fiberと近傍のpartitionをcompact X_d上で保つと半空間germが残る（`ncc:elliptic-face-limitation`）。外側S³への接着だけでこのgermを消すことはできない。これはnormal exteriorのfibrationが存在するという成果と、境界のないglobal baseの未完成を分ける点である。

**相対handles。** Eはouter S³からHopf linkに沿う一つのround 1-handleで得られ、ordinary 1+2 handlesへ分けられる。逆向きは2+3で、唯一の2-handleのbased attaching wordはwhiskerと向きでbにできる。Knot isotopyやframingの自由な選択は主張しない。V≃V₀∪_i(C′×(S¹∨S²))はV₀からの写像を保ち、i=f_e∘q、q:C′→S²、[f_e]=e hである。CofiberはwedgeだがV自体はV₀とのwedgeではない。Curve 2-cell×normal 1-cellがe(τ−1)hの関係を担う。

C′×Eのouter側のminimal relative handle数は(1,2g+1,2g+1,1)、index1–4。H₄=Z、reverseではH₅=Zのため、どちら向きにも全二境界を保つWeinstein cobordismにならない（`tnc:no-weinstein-normal-transition`）。C°のboundary componentsがr≥1ならouter relative Betti数はindex1–3で(1,2g+r,2g+r−1)となる。これは側面∂C°×Eをrelative subspaceに含めない計算であり、puncture後のWeinstein構造の存在は未証明。

**Normal disk module。** `tnc:normal-relative-disk-module` はπ₂(E,Y)=ΛB、∂(fB)=f(1)b+f′(1)a、sphere subgroup=(τ−1)²ΛBを与える。Bは適切に選んだ唯一のreverse 2-handleのcore。任意のb-capではなく、B+(τ−1)²gBは係数がunitでなければΛ-basisにならない。Yの誘導被覆T²×RでFが実際のsolid torusで埋まるためcentral extension classは零。AbeliannessとPeiffer identityが境界kernelの作用を消し、一般relative Hurewiczで必要な作用quotientを正当化する。ν_relは旧sphere generatorとunit ±τⁿだけ異なり、marked orientationとの同一視は別途指定する。

**Curve disksの非可換関係。** `tnc:relative-disk-central-extension` でπ₂(V,∂U)はM=π₂(V)によるπ₁(C′)×Z²のcentral extensionとなり、∏[A_i,B_i]=e h。Curve diskを変更してもcommutator積は変わらず、t作用との整合性がe(τ−1)h=0を再現する。Hyperplane degree ℓでMをZへ押し出すとEuler数eのcircle bundle P_eのπ₁とZ²の積になる。全三次元pairの各curve diskへのt作用、embedded disksとsmooth framingはまだ決定していない。∂UはLagrangianでないので相対群全体へsymplectic areaを定義しない。Sphere subgroupではarea=aℓ、c₁=(5−d)ℓ。

**Clean diskの制約。** Ω=dX∧dYについてμ(a)=0, μ(b)=2（`ncc:normal-maslov-marking`）。Normal C²内のLagrangian diskがregular torusと境界でcleanに交わるなら、境界ma+nbはn=0かつmA=0。逆向き2-handleのword b+kaをこの種のcoreで実現する案は使えない。Elliptic collapseのtransverse diskはsymplecticなので矛盾しない。Curve方向へ動くdiskやnonclean/singular diskの非存在には拡張しない。

**一次資料と位置づけ。** M. Symington, *Four dimensions from two in symplectic topology*, Proc. Sympos. Pure Math. 71 (2003), 153–208, arXiv:math/0210033, §6.1 Lemma 6.3 / Theorem 6.5は標準nodal trade。D. Asimov, *Round handles and non-singular Morse–Smale flows*, Ann. of Math. 102 (1975), 41–54, DOI 10.2307/1970972はround handleの出典（公式書誌を確認、本稿のhandle分解は直接証明）。K. Cieliebak–Y. Eliashberg, *Flexible Weinstein manifolds*, arXiv:1305.1635, §1.1 pp.2–3はisotropic stable manifoldsとindex制限。R. Brown, JLMS (2) 22 (1980), 146–152, DOI 10.1112/jlms/s2-22.1.146, [author PDF](https://www.groupoids.org.uk/pdffiles/whiteheadthm.pdf) §2はrelative crossed module。B. Conrad, [Group cohomology and group extensions](https://math.stanford.edu/~conrad/210BPage/handouts/gpext.pdf), §2はH²によるcentral extension分類。Hatcher *Algebraic Topology*, Props.0.18–0.19、Thms.4.37,4.57はhomotopy pushout、relative Hurewicz、cohomology表現。これらの道具自体と本研究の具体計算を区別し、包括的新規性は未確認。

**次の着手。** Toric outer collarとV₀側を接続し、内部elliptic faceの近傍を変える。C′を有境界曲面へ切る案では、上のH₄障害を除いた後に側面のwhole fibers・action periods・coorientationを定める。Normal diskの二階差分とcurve relationを保ち、global circle actionの障害を避けるgenerator変更領域が必要である。

## 検証
2026-11-21の分割前に相互点検された対象はnormal ballへの明示的補完、elliptic base germ、round handlesとrelative pushout、normal relative crossed module、全補集合のcentral extension、Weinsteinおよびclean normal diskの障害である。原稿の証明と独立点検を行った。既存の補集合位相・short-circle averaging・KE局所解析の証明と点検記録は保持するが、大域接着全体を再監査したものではない。

`python3 verification/verify_focus_relative_models.py --json` を実行し、243 Laurent多項式の境界・整数除法、1,701 deck作用、15 closed-curveと60 punctured-curveのKünneth計算を正確な整数演算で確認した。有限代数検算は相対Hurewiczの仮定・handle geometry・symplectic補完の証明を代替しない。外部査読と形式証明は未実施。

2026-11-21に三次元局所稿190頁・全体稿211頁をLuaLaTeXとlatexmkで組版し、未解決の参照・引用がないことを確認した。新節の代表頁と図はPDFを画像化して表示を点検した。曲面主稿改訂8の130頁は保持した組版記録で、今回再組版していない。既存節の2.27ptのOverfull警告一件は残るが、新節にはない。GitHub反映は組版と別に確認する。