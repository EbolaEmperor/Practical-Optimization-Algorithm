# 实用优化算法 代码库

这是一个使用纯c++编写的实用优化算法库，仍在持续更新。现有内容如下：

- 0.基础运算包
  - 浮点矩阵运算库
  - 高精度整数运算库
  - 高精度分数运算库
  - 高精度分数矩阵运算库
- 1.线搜索
  - Fibonacci分割法
  - Newton单点插值法
  - Armijo-Goldstein不精确线搜索法
  - Wolfe-Powell不精确线搜索法
- 2.最速下降法与牛顿法
  - 最速下降法
  - 牛顿法
  - Gill-Murray稳定牛顿法
- 3.共轭梯度法
  - 共轭梯度法
  - 精确共轭梯度法
  - 预优共轭梯度法
- 4.拟牛顿法
  - 采用简单准则的BFGS校正法
  - 采用Wolfe准则的BFGS校正法（表现最好）
  - 采用Goldestein准则的BFGS校正法（通常不使用）
  - 采用Wolfe准则的DFP校正法
  - 采用Wolfe准则的Broyden族校正法

### 一些说明

目前发现Gill-Murray稳定牛顿法存在bug，暂时不能使用