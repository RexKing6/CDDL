# Zero-Shot-CDDL
研究论文《Zero-Shot Image Classification via Coupled Discriminative Dictionary Learning》已被国际学术会议LSMS-ICSEE2017接受录用，论文代码以FDDL为参考。

main.m 为主函数
count.m 为测试结果计数
各函数作用:
读取数据 ->
参数设置 ->
属性D1,z1初始化并进行第一轮小优化 ->
特征D2,z2初始化并进行第一轮小优化 ->
第二轮大优化 ->
测试

smalldata.mat: 之前为了检测代码是否能运行而提取的小数据
bigdata.mat: 从AwA每一类提取90张图片之后pca降维的数据

FDDL_INIC: 初始化系数z
FDDL_INID: 初始化字典D
FDDL_UpdateDi: 更新字典
Initround1: 初始化并进行第一轮优化
IPM_SC: 根据数据和字典计算出系数z
Round1_Class_Energy: 第一轮内部每一类最小值计算
Round1_FDL_Energy: 第一轮全部类最小值计算
Round1_Gradient_Comp: 计算论文算法中的倒三角Q
Round1_SpaCoef: 第一轮更新系数
Round2: 第二轮大优化主函数
Round2_Class_Energy: 第二轮内部每一类最小值计算
Round2_FDL_Energy: 第二轮全部类最小值计算
Round2_Gradient_Comp: 计算论文算法中的倒三角Q
Round2_SpaCoef: 第二轮更新系数
soft: 论文算法自带的soft函数
