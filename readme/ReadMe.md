[TOC]

# 程序简介

本程序采用python3编写，具体用到的第三方库有numpy，pyyaml，scipy。如有问题，及时联系J.Pei(J.Pei@foxmail.com)。

# 程序准备条件

运行程序前需要明确以下信息:

1. 纵波声速，横波声速 (奥林巴斯超声回声测试仪可测)
2. 样品密度(阿基米德排水法或者**Jade**分析XRD可知)
3. 单胞中原子个数和原子体积(**ICSD**数据库或者**Jade**分析XRD可知)
4. 基体的晶格热导率与温度的关系，并写入到”**xxx.dat**“文件中。”**xxx.dat**“文件每行前有 *#* 表示注释，程序不会读取此行数据。温度与热导率之间需至少一个空格。”**xxx.dat**“文件格式如下:

```dat
##
303.2	0.840340632511578
324.5	0.806283215996452
349.2	0.788495783531084
375.6	0.783078527714649
400	0.792257053666418
423.1	0.817909863466827
445.6	0.85217890658237
471	0.907103882898424
496.4	0.977540942490175
522.6	1.06148964846849
548.4	1.16173446302351
573.7	1.24443211859626
```

> ”**xxx.dat**“请用"**notepad++**"或者“**gedit**"等软件写的，一定不能用**Windows**下的”**记事本**“，会报错。

5. 固溶体的化学表达式。

# 程序原理公式

## 平均声速[^07-04]

$$
v_{average} = \left(\frac{1}{3}\left[\frac{1}{v_{l}^{3}}+\frac{2}{v_{s}^{3}}\right]\right)^{-\frac{1}{3}}
$$

其中，$v_{average}$是平均声速，$v_{l}$表示纵波声速,$v_{s}$表示横波声速。

[^07-04]: K. Kurosaki, A. Kosuga, H. Muta, M. Uno, S. Yamanaka, Ag9TlTe5: A high-performance thermoelectric bulk material with extremely low thermal conductivity, Appl. Phys. Lett. 87 (2005) 3–6. doi:10.1063/1.2009828.



德拜温度[^07-04]

$$
\theta_{D}=\frac{h}{k_{B}}\left[\frac{3N}{4\pi V}\right]^{1/3}v_{average}
$$

其中，$h$为普朗克常量，$k_{B}$为玻尔兹曼常数，$N$为单胞中的原子数目，$V$为单胞体积。

## 德拜温度

- 对于各向同性的材料来说[^04-04]

$$
\theta_D=\left( \frac {h}{k_{B}} \right) \left[ \frac{9N}{4\pi V(v_{L}^{-3}+2v_{S}^{-3})}\right]^{1/3}
$$

其中$\theta_D$代表德拜温度，$h$代表普朗克常量，$k_B$代表玻尔兹曼常数，$N$代表单胞中的原子数目，$V$代表晶胞体积，$v_{L}$代表纵向声速，$v_{S}$代表剪切声速

[^04-04]: A. Kosuga, M. Uno, K. Kurosaki, and S. Yamanaka, “Thermoelectric properties of Ag1−xPb18SbTe20 (x = 0, 0.1, 0.3),” J. Alloys Compd., vol. 387, no. 1–2, pp. 52–55, 2005.

- 德拜温度的实验反推[^04-11]

通过PPMS测试，测得低温下1.8K-300K范围内的低温热容$C_{p}$，通过$C_{p}/T$与$T^{2}$的线性拟合关系，可以推出德拜温度，这里有两种模型，一种是德拜模型，一种是德拜-爱因斯坦模型。



## 热导率的非晶极限

- Cahill公式

Cahill公式适用的是非晶材料热导率的高温极限(最小热导率不考虑电子热导，仅从晶格热导的角度出发)。

- 简化公式[^04-09]

$$
\kappa_{min}=\frac{1}{2}\left(\frac{\pi}{6}\right)^{1/3}k_{B}V^{-2/3}(2v_{t}+v_{l})
$$



其中，$\kappa_{min}$代表热导率的非晶极限, $k_{B}$代表玻尔兹曼常数，$V$代表平均原子体积，$v_{l}$代表纵波声速，$v_{t}$代表横波声速。

[^04-09]: A. F. May and G. J. Snyder, “Introduction to Modeling Thermoelectric Transport at High Temperatures,” Thermoelectr. Its Energy Harvest., p. 11, 2012.

> 注:谭刚健[^04-14]等人用的公式与上上不同.具体公式如下:

$$
\kappa_{min} = \frac{\pi}{4}k_{B}V^{-2/3}\nu
$$

其中，$\nu$表示平均声速。

[^04-14]: Tan G, Shi F, Doak JW, Sun H, Zhao L-D, Wang P, et al. Extraordinary role of Hg in enhancing the thermoelectric performance of p-type SnTe. Energy Environ Sci 2015;8:267–77. doi:10.1039/C4EE01463D.

## 杨氏模量$E$

- 对于各向同性的材料来说：[^07-01]


$$
E=\frac{G(3v_{L}^{2}-4v_{S}^{2})}{v_{L}^{2}-v_{s}^{2}}
$$

其中$E$代表杨氏模量,$G$代表剪切模量，$v_{L}$代表纵向声速，$v_{S}$代表剪切声速

[^07-01]: A. Kosuga, M. Uno, K. Kurosaki, and S. Yamanaka, “Thermoelectric properties of Ag1−xPb18SbTe20 (x = 0, 0.1, 0.3),” J. Alloys Compd., vol. 387, no. 1–2, pp. 52–55, 2005.

- 基于弹性常数的计算[^07-02]

$$
 E=\frac{C_{44}(3C_{11}-4C_{44})}{C_{11}-C{44}}
$$

其中，$E$代表杨氏模量，$C_{11}，C_{44}$代表不同方向的弹性常数 。

[^07-02]: F. Ren, E. D. Case, J. E. Ni, E. J. Timm, E. Lara-Curzio, R. M. Trejo, C. H. Lin, and M. G. Kanatzidis, “Temperature-dependent elastic moduli of lead telluride-based thermoelectric materials,” Philos. Mag., vol. 89, no. 2, pp. 143–167, 2009.

## 剪切模量$G$

- 对于各向同性的材料来说：[^07-01]

$$
G=\rho v_{S}^{2}
$$

其中，$G$代表剪切模量，$\rho$代表样品的密度，$v_{S}$代表剪切声速

## 块体模量$B$计算[^07-03]

$$
B=\rho\left(v_{L}^{2}-\frac{4}{3}v_{T}^{2}\right)
$$

通过声速测试得到纵波声速$v_{L}$和横波声速$v_{T}$，再乘以密度$\rho$即可得到块体模量$B$。

[^07-03]: G. Tan et al., “High Thermoelectric Performance in SnTe-AgSbTe2Alloys from Lattice Softening, Giant Phonon-Vacancy Scattering, and Valence Band Convergence,” ACS Energy Lett., vol. 3, no. 3, pp. 705–712, 2018.

## 体积模量K计算

$$
K=\frac{E}{3(1-2\nu)}
$$

其中，$E$为杨氏弹性模量，$\nu$为泊松比，本公式选择材科基。

## 泊松比

- 对于各向同性的材料来说
  - 基于弹性常数的计算[^07-02]

$$
\nu =\frac{C_{11}-2C_{44}}{2(C_{11}-C_{44})}
$$

其中，$\nu$代表泊松比，$C_{11}，C_{44}$代表不同方向的弹性常数
$$
\nu = \frac{\left(\frac{v_{L}}{v_{T}}\right)^{2}-2}{2 \left[ \left(\frac{v_{L}}{v_{T}}\right)^{2}-1 \right]}
$$

## 格林艾森常数

$$
\gamma=\frac{3}{2}\left(\frac{1+\nu}{2-3\nu}\right)
$$



## 应力场相关的适配参数$\varepsilon$

$$
\varepsilon=\frac{2}{9}\cdot\left(\frac{6.4\times\gamma\times(1+\nu)}{1-\nu}\right)^{2}
$$



## deby-callaway模型[^07-05][^07-06][^07-07][^07-08]

$$
\frac{\kappa_{L}}{\kappa_{L,p}}=\frac{tan^{-1} u}{u}
$$

$$
u^{2}=\frac{\pi^{2}\theta_{D}\Omega}{hv_{average}^{2}}\kappa_{L,p}\Gamma
$$

$$
\Gamma = \Gamma_{S}+\Gamma_{M}
$$

$$
\Gamma_{M} = \frac{\sum_{i=1}^{n}c_{i}\left(\frac{\overline{M}_{i}}{\widehat{M}}\right)^{2}f_{i,1}f_{i,2}\left(\frac{M_{i,1}-M_{i,2}}{\overline{M}_{i}}\right)^{2}}{\sum_{i=1}^{n}c_{i}}
$$

$$
\Gamma_{S} = \frac{\sum_{i=1}^{n}c_{i}\left(\frac{\overline{M}_{i}}{\widehat{M}}\right)^{2}f_{i,1}f_{i,2}\varepsilon_{i}\left(\frac{r_{i,1}-r_{i,2}}{\overline{r}_{i}}\right)^{2}}{\sum_{i=1}^{n}c_{i}}
$$

$$
\overline{M}_{i}=\sum_{k}f_{i,k}M_{i,k}
$$

$$
\overline{r}_{i}=\sum_{k}f_{i,k}r_{i,k}
$$

$$
\widehat{M}=\frac{\sum_{i=1}^{n}c_{i}\overline{M}_{i}}{\sum_{i=1}^{n}c_{i}}
$$

其中，$\kappa_{L}$为固溶体中计算的晶格热导率,$\kappa_{L,p}$为基体的晶格热导率(可由SPB模型计算得),$u$为无序度因子，$\theta_{D}$为德拜温度,$\Omega$为平均原子体积，$h$为普朗克常量,$v_{average}$为平均声速,$\Gamma$为散射参数，可以由质量波动$\Gamma_{M}$和应力波动$\Gamma_{S}$求得。$n$为主晶格中不同原子位置个数，以(Bi,Sb)2Te3 或者GeTe为例，n=2。$c_{i}$为原子占位的简并度，以GeTe为例，$c_{1}=c_{2}=1$;以(Bi,Sb)2Te3为例，$c_{1}=2;\; c_{2}=3$。$\widehat{M}$为该化合物的平均原子质量，$\overline{M}_{i}$和$\overline{r}_{i}$为第ith子晶格的平均原子质量和原子半径，$f_{i,k}$为第ith子晶格中第kth个原子的占比分数，以$Fe_{0.3}(Bi,Sb)_{1.7}Te_{3}$为例，Fe的$f_{i,k}=0.15$, (Bi,Sb)的$f_{i,k}=0.85$。

[^07-05]: 肖钰博士论文
[^07-06]: 姚瑶博士论文
[^07-07]: Zheng, Z., Su, X., Deng, R., Stoumpos, C. C., Xie, H., Liu, W., … Tang, X. (2018). Rhombohedral to Cubic Conversion of GeTe via MnTe alloying Leads to Ultralow Thermal Conductivity, Electronic Band Convergence and High Thermoelectric Performance. *Journal of the American Chemical Society*, jacs.7b13611. https://doi.org/10.1021/jacs.7b13611
[^07-08]: Pan, Y., & Li, J.-F. (2016). Thermoelectric performance enhancement in n-type Bi2(TeSe)3 alloys owing to nanoscale inhomogeneity combined with a spark plasma-textured microstructure. *NPG Asia Materials*, *8*(6), e275. https://doi.org/10.1038/am.2016.67

# 参考文献

