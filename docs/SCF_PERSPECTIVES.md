# SCF 的两个视角

## 一、迭代视角:解一个不动点方程

自洽场的难处是一个循环依赖。Fock 矩阵由密度决定,而密度又由 Fock 矩阵的本征轨道填出来:

$$
\mathbf F[\mathbf P]\,\mathbf C=\mathbf S\,\mathbf C\,\boldsymbol\varepsilon,
\qquad
\mathbf P=\mathbf C_{\rm occ}\mathbf C_{\rm occ}^{T}.
$$

把"由 $\mathbf P$ 造 $\mathbf F$、对角化、再填出新 $\mathbf P$"这一整套记作映射 $\mathbf P\mapsto g(\mathbf P)$,SCF 就是求

$$
\mathbf P=g(\mathbf P).
$$

最朴素的做法是直接代入 $\mathbf P_{n+1}=g(\mathbf P_n)$,但这个迭代常常振荡甚至发散,所以要做加速。收敛时 $\mathbf F$ 与 $\mathbf P$ 在 $\mathbf S$ 度量下对易,于是定义误差矩阵

$$
\mathbf e=\mathbf F\mathbf P\mathbf S-\mathbf S\mathbf P\mathbf F .
$$

DIIS 记住最近 $m$ 步的 $(\mathbf F_i,\mathbf e_i)$,在 $\sum_i c_i=1$ 的约束下最小化组合误差的范数,再用同一组系数拼出下一个 Fock:

$$
\min_{\{c_i\}}\Big\|\sum_i c_i\,\mathbf e_i\Big\|^2,
\qquad
\mathbf F_{n+1}^{\rm extr}=\sum_i c_i\,\mathbf F_i .
$$

远离收敛时,DIIS 的线性假设不可靠,ADIIS 改为最小化能量的一个近似模型来选 $c_i$。这一视角的特点是每圈只需要一次 Fock 构建,但它只关心"怎么让不动点迭代收敛",对能量面本身的形状一无所知,所以在近简并、金属和振荡体系上会卡住。

## 二、轨道旋转视角:在能量面上找最低点

换一个问法:能量是轨道的函数,SCF 就是求它的极小值。保持轨道正交的最自然的参数化,是把占据轨道向虚轨道混合。取反对称的旋转生成元 $\boldsymbol\kappa$,它只有占据–虚的块 $\kappa_{ai}$(共 $n_{\rm occ}n_{\rm vir}$ 个自由度),则

$$
\mathbf C(\boldsymbol\kappa)=\mathbf C_0\,e^{\boldsymbol\kappa},
\qquad
\boldsymbol\kappa=\begin{pmatrix}0&-\boldsymbol\kappa^{T}\\ \boldsymbol\kappa&0\end{pmatrix}.
$$

这样 SCF 变成关于 $\boldsymbol\kappa$ 的无约束极小化。对能量在 $\boldsymbol\kappa=0$ 处展开到二阶:

$$
E(\boldsymbol\kappa)\approx E_0+\mathbf g^{T}\boldsymbol\kappa+\tfrac12\,\boldsymbol\kappa^{T}\mathbf H\,\boldsymbol\kappa .
$$

一阶项就是 Fock 矩阵的占据–虚块,收敛条件 $\mathbf g=0$ 正是 Brillouin 条件:

$$
g_{ai}=4\,F_{ai}\ (\text{闭壳层}),\qquad g_{ai}=0 .
$$

二阶项为轨道 Hessian,它由轨道能差和一个响应项组成,与 TDDFT 的 $(\mathbf A+\mathbf B)$ 是同一个算符:

$$
H_{ai,bj}=\delta_{ab}\delta_{ij}(\varepsilon_a-\varepsilon_i)\cdot 4
\;+\;\big[(\mathbf A+\mathbf B)-\text{对角部分}\big]_{ai,bj},
$$

其中响应部分包含 Coulomb、$-c_x$ 倍的精确交换,以及 $f_{\rm xc}$ 核。

在这个框架下,有两种取步策略。

### 策略 A:只用梯度(一阶与拟牛顿)

沿下降方向 $\mathbf d=-\mathbf M^{-1}\mathbf g$ 走一步,步长 $\alpha$ 由线搜索决定:

$$
\boldsymbol\kappa\leftarrow\alpha\,\mathbf d .
$$

最简单的 $\mathbf M$ 取 Hessian 的对角部分 $4(\varepsilon_a-\varepsilon_i)$。SOSCF 一类方法再用最近几步的 $\mathbf s_k=\Delta\boldsymbol\kappa,\ \mathbf y_k=\Delta\mathbf g$ 通过 L-BFGS 递推去补非对角信息,满足割线条件

$$
\mathbf H_{k+1}\,\mathbf s_k=\mathbf y_k .
$$

每步只要一次 Fock 构建,能量沿搜索方向单调下降,但后期是线性收敛。

### 策略 B:用 Hessian(牛顿与信赖域)

直接解二次模型的驻点方程

$$
\mathbf H\,\boldsymbol\kappa=-\mathbf g .
$$

$\mathbf H$ 的维数是 $n_{\rm occ}n_{\rm vir}$ 的平方,存不下,所以用共轭梯度(CG)迭代,只需要 Hessian 乘向量

$$
\mathbf H\mathbf v=\big[4(\varepsilon_a-\varepsilon_i)v_{ai}\big]+\mathbf J[\mathbf v]-c_x\mathbf K[\mathbf v]+f_{\rm xc}\!\cdot\!\mathbf v .
$$

这个乘积也可以用有限差分得到,只多做一次 Fock 构建:

$$
\mathbf H\mathbf v\approx\frac{\mathbf g(\boldsymbol\kappa=h\mathbf v)-\mathbf g(0)}{h}.
$$

信赖域版本把步长限制在半径 $\Delta$ 内,保证即使在 Hessian 有负本征值(SCF 解不稳定)时也有定义:

$$
\min_{\|\boldsymbol\kappa\|\le\Delta}\ \mathbf g^{T}\boldsymbol\kappa+\tfrac12\boldsymbol\kappa^{T}\mathbf H\boldsymbol\kappa .
$$

牛顿步是二次收敛的,通常 5 到 8 个宏观步,但每个宏观步需要若干次 $\mathbf H\mathbf v$。

## 三、成本怎么算

两个视角的差别,不在圈数,而在总的 Fock(或 $\mathbf K$)构建次数:

$$
T\;\approx\;N_{\rm Fock}\times t_{\rm Fock},\qquad
N_{\rm Fock}^{\rm DIIS}=N_{\rm cyc},\quad
N_{\rm Fock}^{\rm Newton}=N_{\rm macro}\,(1+N_{\rm CG}).
$$

在 Direwolf 里,二茂铁 B3LYP/def2-SVP RIJCOSX 的实测是:DIIS 18 圈、约 20 次 Fock 构建、29.2 s;Newton-CG 13 个宏观步,共 27 到 41 次构建、48.8 s。圈数少了,但每个 $\mathbf H\mathbf v$ 近似是一次完整的 $\mathbf K$ 构建,所以总时间更长。

因此实践中的分工是:大多数体系用 DIIS 加 ADIIS,又快又稳;DIIS 振荡或不收敛的难体系(过渡金属开壳层、近简并、拉伸键)才交给轨道旋转的方法,策略 A 求稳,策略 B 求最后几位有效数字的二次收敛。
