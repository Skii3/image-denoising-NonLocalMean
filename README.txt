主函数：svd_lra.m
svd实现函数：func_svd_lra4.m

主要算法思路介绍：

1、将二维图像数据转化为patch数据：从图像中按步长取m*m大小，拉成一个向量，总共得到M个向量；
2、对M个patch进行聚类，按照一定步长取中心patch，在该中心patch周围一定大小的矩形框内挑选与之欧氏距离最近的L个patch，构成patch group；
3、对每个patch group进行svd分解，根据噪声大小取前r个特征值，重构得到新的patch矩阵；
4、对所有新的patch group进行聚合，处理重复变换的patch以及像素，恢复出整张图；
5、得到第一次去噪后的图，按一定比例与原始噪声图相加，再进行1~3的处理。
