# 实用优化算法 代码库

这是一个使用纯c++编写的实用优化算法库，仍在持续更新。现有内容如下：

- 0.基础运算包
  - 浮点矩阵运算库 `matrix.h`
  - 高精度整数运算库 `bigint.h`
  - 高精度分数运算库 `fraction.h`
  - 高精度分数矩阵运算库 `matrix_fraction.h`
  - 数值方法求导运算库 `derivation.h`
- 1.线搜索 (`linear_search.h`)
  - Fibonacci分割法
  - Newton单点插值法
  - Armijo-Goldstein不精确线搜索法
  - Wolfe-Powell不精确线搜索法
- 2.最速下降法与牛顿法 (`steepest_descent.h`, `newton.h`)
  - 最速下降法
  - 牛顿法
  - Gill-Murray稳定牛顿法
- 3.共轭梯度法 (`conjugate_gradient.h`, `conjugate_gradient_accurate.h`)
  - 共轭梯度法
  - 重启动FR非线性共轭梯度法
  - 精确共轭梯度法
  - 预优共轭梯度法
- 4.拟牛顿法 (`quassi_newton.h`)
  - 采用简单准则的BFGS校正法
  - 采用Wolfe准则的BFGS校正法（表现最好）
  - 采用Goldestein准则的BFGS校正法（通常不使用）
  - 采用Wolfe准则的DFP校正法
  - 采用Wolfe准则的Broyden族校正法
- 5.信赖域法 (`trust_region.h`)
  - 采用Dogleg方法的信赖域法
- 6.非线性最小二乘问题 (`least_square.h`)
  - Levenberg-Marquardt方法
- 7.线性规划问题 (`linprog.h`)
  - 单纯形方法
- 8.凸二次规划问题 (`quadric_programming.h`)
  - 求解带等式约束的凸二次规划问题的拉格朗日乘子法
  - 求解一般凸二次规划的有效集方法
- 9.一般约束最优化问题 (`general_constraint.h`)
  - 求解一般约束优化问题的乘子法（PHR算法）

有一些方法有原始版本和对应的`gradfree`版本，后者利用数值方法自动求导，当待优化函数由人工求导过于复杂时可以使用。但对于方便求导的函数，还是建议使用原始版本，因为虽然数值方法求导对计算精度影响不大，但速度较慢，影响效率。

## 编译方式

运行下面的命令编译已有 example：
```bash
g++ example/file.cpp -std=c++20 -I. -o test
```
把 `file.cpp` 改成你要编译的 example 的文件名，你会在当前目录下得到一个叫 `test` 的可执行文件。

## 使用模式

提供了两个编译选项：`DEBUG`和`SILENCE`，使用宏定义可以选择模式. 其中通常模式下只输出最优化结果，`DEBUG`模式提供较为详细的中间步骤信息，`SILENCE`模式将不输出任何信息与结果.

## 测试函数集

在func文件夹下有以下几个常用测试函数，直接引用函数名对应的`.h`文件即可，每个文件都包含函数`f(x)`、梯度`grad(x)`、Hessian矩阵`hessian(x)`、初始迭代位置`initial()`

### 高维Rosenbrock函数

$$
f(x)=\sum_{i=1}^{n-1} [100(x_{i+1}-x_i^2)^2+(1-x_i)^2]
$$

### Powell奇异函数

$$
f(x)=(x_1+10x_2)^2+5(x_3-x_4)^2+(x_2-2x_3)^4+10(x_1-x_4)
$$

### Beale函数

$$
f(x) = (1.5-x_1+x_1x_2)^2 + (2.25-x_1+x_1x_2^2)^2 + (2.625-x_1+x_1x_2^3)^2
$$