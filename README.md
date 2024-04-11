# school_work-matrix

**研一课程 矩阵分析与应用 期末作业**

### 程序介绍

出于十进制小数一般无法用二进制精确表示的考虑，将矩阵的元素设为新的自定义类型 Value，它由两个 double 组成，分别表示分子和分母，并可以像浮点数一样使用各种常见的运算符。这种方法保证了在对可逆矩阵进行 LU 分解时，由于每一步只用到了乘除法，因此不会有精度损失，但在进行其他分解时，由于用到了根号运算，仍可能产生累计误差并导致结果错误。LU 分解，QR 分解，Householder 约简，Givens 约简都能得到三角矩阵，可用于快速求解线性方程组，但仅限于方程组有唯一特解或解不存在时求它的最小二乘解。求矩阵行列式是通过 Givens 约简实现的，在 PA=T 中，P 是旋转矩阵的乘积，且为正交阵，因此 det(P)=1，det(T)为对角元素的乘积，det(A)=det(T)。

### 程序运行说明和实例

在对矩阵进行操作前，需要先输入矩阵，否则当矩阵为空时，各项操作都无法进行。输入矩阵时，首先输入维度 m 和 n，然后根据提示依次输入每个元素。可在一行输入多个元素，每个元素的形式是“a/b”或“a”，a 和 b 均为浮点数。下图所示为输入一个 3x3 矩阵的过程。

![](images/1-image.png)

输入完矩阵后，若A不为空，则会显示它的“状态”，会计算矩阵的秩，并显示能对矩阵进行的分解操作。线性方程组的右端b的维度会随矩阵A的变化而变化，即b的行数始终与A相同，列数始终为1，输入b时无需再设置维度。下图所示为输入b的过程。

![](images/2-image.png)

若矩阵A的行数改变会重新设置b的元素全为0. 可依据提示对A进行各项分解，若分解能得到三角矩阵，则还会自动给出线性方程组的解。下图所示为LU分解的结果。

![](images/3-image.png)

类似的，输入4，进行QR分解。由于计算误差的原因，结果有些许不同。

![](images/4-image.png)

输入5，进行Householder约简。

![](images/5-image.png)

输入6，进行Givens约简。结果应与Householder约简的结果一致。

![](images/6-image.png)

输入7，进行URV分解，此处计算产生了溢出，仍需要进行改进。

![](images/7-image.png)

输入8，计算矩阵A的行列式，当A不可逆时，行列式为0，A可逆时使用Givens约简计算行列式。根据命令前给出的对A的描述可知道A是否有行列式。

![](images/8-image.png)

当输入0时退出程序，输入字母也会被当作0，会退出程序。

![](images\9-image.png)
