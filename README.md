# CDDL
《Zero-Shot Image Classification via Coupled Discriminative Dictionary Learning》LSMS-ICSEE2017。论文代码很大程度参考[FDDL](http://www4.comp.polyu.edu.hk/~cslzhang/papers.htm)和[SSE](https://zimingzhang.wordpress.com/source-code/)。

## Dateset

* AwA(https://cvml.ist.ac.at/AwA/)
* aP&Y(http://vision.cs.uiuc.edu/attributes/)

## Requirements

* Matlab==2015a

## Usage

1. Download the feature mat from https://drive.google.com/drive/folders/1VzE84JJdnl45Iqh-KBPod0ZE5Zt5b7se, and then put it in `./`
2. run `main.m`

## Description

1. 读取数据
2. 参数设置
3. 属性D1,z1初始化并进行第一轮小优化
4. 特征D2,z2初始化并进行第一轮小优化
5. 第二轮大优化
6. 测试

## Functions

* `main.m`: 主函数
* `count.m`: 测试结果计数
* `FDDL_INIC.m`: 初始化系数z
* `FDDL_INID.m`: 初始化字典D
* `FDDL_UpdateDi.m`: 更新字典
* `Initround1.m`: 初始化并进行第一轮优化
* `IPM_SC.m`: 根据数据和字典计算出系数z
* `Round1_Class_Energy.m`: 第一轮内部每一类最小值计算
* `Round1_FDL_Energy.m`: 第一轮全部类最小值计算
* `Round1_Gradient_Comp.m`: 计算论文算法中的$\nabla Q$
* `Round1_SpaCoef.m`: 第一轮更新系数
* `Round2.m`: 第二轮大优化主函数
* `Round2_Class_Energy.m`: 第二轮内部每一类最小值计算
* `Round2_FDL_Energy.m`: 第二轮全部类最小值计算
* `Round2_Gradient_Comp.m`: 计算论文算法中的$\nabla Q$
* `Round2_SpaCoef.m`: 第二轮更新系数
* `soft.m`: 论文算法自带的soft函数

## References 

* Farhadi A, Endres I, Hoiem D, et al. Describing objects by their attributes[C]// Computer
  Vision and Pattern Recognition, 2009. CVPR 2009. IEEE Conference on. IEEE, 2009:1778-1785.
* Lampert C H, Nickisch H, Harmeling S. Attribute-Based Classification for Zero-Shot Visual
  Object Categorization[J]. IEEE Transactions on Pattern Analysis & Machine Intelligence,
  2014, 36(3):453-65. 
* Yang M, Zhang L, Feng X, et al. Fisher Discrimination Dictionary Learning for sparse representation[C]// International Conference on Computer Vision. IEEE, 2011:543-550. 
* Zhang Z, Saligrama V. Zero-Shot Learning via Semantic Similarity Embedding[C]// IEEE
  International Conference on Computer Vision. IEEE, 2015:4166-4174. 
