---
title: R 中的地理空间插值
date: 2023-02-18 12:34:12 +0800
categories: [R, 空间插值]
tags: [空间插值, R, 空间分析]     # TAG names should always be lowercase
---


克里格(Kriging)是一组地质统计学技术，用于从附近位置的已知观测值中插值非采样位置的随机场值(例如，土壤变量作为地理位置的函数)。克里格法背后的主要统计假设是平稳性，这意味着统计属性(如均值和方差)不依赖于确切的空间位置，因此一个变量在一个位置的均值和方差等于另一个位置的均值和方差。克里格的基本思想是通过计算该点附近函数已知值的加权平均来预测给定点上的值。与其他确定性插值方法(如反向距离加权(IDW)和样条)不同，克里格基于自相关(即测量点之间的统计关系)来插值空间场中的值。Kriging能够生成不确定的预测曲面。虽然平稳性(恒定的均值和方差)和各向同性(所有方向的均匀性)是克里格提供最佳线性无偏预测的两个主要假设，但是，这些假设对于各种形式和方法的克里格有灵活性。

# 克里格法的步骤:

- 搜索半径内样本的平均样本对样本变异性的计算；
- 选择位于搜索半径内的最近样本；
- 建立kriging矩阵，包括建立包含每个邻域样本值的预期样本变异性的半方差矩阵，以及建立包含每个最接近邻域样本值与块之间的平均变异性的矩阵；
- 克里格系数矩阵的建立; 以及将克里格系数乘以它们各自的样本值，以提供克里格估计 (KE)。克里格方差 (KV) 是根据权重系数的乘积和它们各自的样本块方差计算得出的。一个额外的常数，滞后范围乘数被添加以最小化方差。

## 普通克里金法

```R
# Variogram
v<-variogram(SOC.bc~ 1, data = train, cloud=F)
# Intial parameter set by eye esitmation
m<-vgm(1.5,"Exp",40000,0.5)
# least square fit
m.f<-fit.variogram(v, m)
m.f
```
### 绘制变量图和拟合模型

```R
#### Plot varigram and fitted model:
plot(v, pl=F, 
     model=m.f,
     col="black", 
     cex=0.9, 
     lwd=0.5,
     lty=1,
     pch=19,
     main="Variogram and Fitted Model\n Box-Cox Transformed SOC",
     xlab="Distance (m)",
     ylab="Semivariance")
```

![xcxc](../../assets/img/test/3_1_9.png)

