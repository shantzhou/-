# Harris特征点检测器-兴趣点检测

### Harris角点检测算法原理

##### 1.算法思想

算法的核心是利用局部窗口在图像上进行移动，判断灰度是否发生较大的变化。如果窗口内的灰度值（在梯度图上）都有较大的变化，那么这个窗口所在区域就存在角点。

这样就可以将 Harris 角点检测算法分为以下三步：

- 当窗口（局部区域）同时向 x （水平）和 y（垂直） 两个方向移动时，计算窗口内部的像素值变化量 $E(x,y)$ ；
- 对于每个窗口，都计算其对应的一个角点响应函数 $R$；
- 然后对该函数进行阈值处理，如果 $R > threshold$，表示该窗口对应一个角点特征。

### 2.Step1  建立数学模型 

**第一步是通过建立数学模型，确定哪些窗口会引起较大的灰度值变化。** 让一个窗口的中心位于灰度图像的一个位置$(x,y)$，这个位置的像素灰度值为$I(x,y)$ ，如果这个窗口分别向 $x$ 和 $y$ 方向移动一个小的位移$u$和$v$，到一个新的位置 $(x+u,y+v)$ ，这个位置的像素灰度值就是$I(x+u,y+v)$ 。

$|I(x+u,y+v)-I(x,y)|$就是窗口移动引起的灰度值的变化值。

设$w(x,y)$为位置$(x,y)$处的窗口函数，表示窗口内各像素的权重，最简单的就是把窗口内所有像素的权重都设为1，即一个均值滤波核。

当然，也可以把 $w(x,y)$设定为以窗口中心为原点的高斯分布，即一个高斯核。如果窗口中心点像素是角点，那么窗口移动前后，中心点的灰度值变化非常强烈，所以该点权重系数应该设大一点，表示该点对灰度变化的贡献较大；而离窗口中心（角点）较远的点，这些点的灰度变化比较小，于是将权重系数设小一点，表示该点对灰度变化的贡献较小。

则窗口在各个方向上移动 $(u,v)$所造成的像素灰度值的变化量公式如下：

![68747470733a2f2f696d67636f6e766572742e6373646e696d672e636e2f6148523063484d364c7939706257466e5a584d794d4445314c6d4e75596d78765a334d75593239744c324a73623263764e4455784e6a59774c7a49774d5459774e4338304e5445324e6a41744d6a41784e6a41304d6](D:\github\图像处理下\68747470733a2f2f696d67636f6e766572742e6373646e696d672e636e2f6148523063484d364c7939706257466e5a584d794d4445314c6d4e75596d78765a334d75593239744c324a73623263764e4455784e6a59774c7a49774d5459774e4338304e5445324e6a41744d6a41784e6a41304d6.jpg)

若窗口内是一个角点，则$E(u,v)$的计算结果将会很大。

为了提高计算效率，对上述公式进行简化，利用泰勒级数展开来得到这个公式的近似形式：

对于二维的泰勒展开式公式为：
$T(x,y)=f(u,v)+(x-u)f_x(u,v)+(y-v)f_y(u,v)+....$

则$I(x+u,y+v)$ 为：
$I(x+u,y+v)=I(x,y)+uI_x+vI_y$

其中$I_x$和$I_y$是$I$的微分（偏导），在图像中就是求$x$ 和 $y$ 方向的**梯度图**：

$I_x=\frac{\partial I(x,y)}{\partial x}$

$I_y=\frac{\partial I(x,y)}{\partial y}$

将$I(x+u,y+v)=I(x,y)+uI_x+vI_y$代入$E(u，v)$可得：

![68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333830383433342e706e67](D:\github\图像处理下\68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333830383433342e706e67.png)

提出 u 和 v ，得到最终的近似形式：

![68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333233333536342e706e67](D:\github\图像处理下\68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333233333536342e706e67.png)

其中矩阵M为：

![68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333235383134352e706e67](D:\github\图像处理下\68747470733a2f2f696d672d626c6f672e6373646e696d672e636e2f32303230303631303132333235383134352e706e67.png)

最后是把实对称矩阵对角化处理后的结果，可以把R看成旋转因子，其不影响两个正交方向的变化分量。

经对角化处理后，将两个正交方向的变化分量提取出来，就是 λ1 和 λ2（特征值）。 这里利用了**线性代数中的实对称矩阵对角化**的相关知识，有兴趣的同学可以进一步查阅相关资料。

### 3.step 2  角点响应函数R

现在我们已经得到 $E(u,v)$的最终形式，别忘了我们的目的是要找到会引起较大的灰度值变化的那些窗口。

灰度值变化的大小则取决于矩阵M，M为梯度的协方差矩阵。在实际应用中为了能够应用更好的编程，所以定义了角点响应函数R，通过判定R大小来判断像素是否为角点。

计算每个窗口对应的得分（角点响应函数R定义）：

![3.1](D:\github\图像处理下\3.1.png)

其中 $det(M)=\lambda_1\lambda_2$是矩阵的行列式， $trace(M)=\lambda_1+\lambda_2$ 是矩阵的迹。

$λ1$ 和 $λ2$ 是矩阵$M$的特征值， $k$是一个经验常数，在范围 (0.04, 0.06) 之间。

$R$的值取决于$M$的特征值，对于角点$|R|$很大，平坦的区域$|R|$很小，边缘的$R$为负值。

### 4.step3 角点判定 

根据 R 的值，将这个窗口所在的区域划分为平面、边缘或角点。为了得到最优的角点，我们还可以使用非极大值抑制。

注意：Harris 检测器具有旋转不变性，但不具有尺度不变性，也就是说尺度变化可能会导致角点变为边缘。想要尺度不变特性的话，可以关注SIFT特征。

因为特征值 λ1 和 λ2 决定了 R 的值，所以我们可以用特征值来决定一个窗口是平面、边缘还是角点：

- 平面:：该窗口在平坦区域上滑动，窗口内的灰度值基本不会发生变化，所以 $|R|$ 值非常小，在水平和竖直方向的变化量均较小，即 $I_x$和 $I_y$都较小，那么 λ1 和 λ2 都较小；
- 边缘：$|R|$值为负数，仅在水平或竖直方向有较大的变化量，即 $I_x$和 $I_y$只有一个较大，也就是 λ1>>λ2 或 λ2>>λ1；
- 角点：[公式] 值很大，在水平、竖直两个方向上变化均较大的点，即 $I_x$和 $I_y$ 都较大，也就是 λ1 和 λ2 都很大。



如下图所示：

![4.1](D:\github\图像处理下\4.1.jpg)



Harris 角点检测的结果是带有这些分数 R 的灰度图像，设定一个阈值，分数大于这个阈值的像素就对应角点。



## 代码实现

Opencv中有提供Harris角点检测的函数，

函数原型：`cv2.cornerHarris(src, blockSize, ksize, k[, dst[, borderType]])`

对于每一个像素 (x,y)，在 (blockSize x blockSize) 邻域内，计算梯度图的协方差矩阵 $M(x,y)$，然后通过上面第二步中的角点响应函数得到结果图。图像中的角点可以为该结果图的局部最大值。

即可以得到输出图中的局部最大值，这些值就对应图像中的角点。

参数解释：

- src - 输入灰度图像，float32类型
- blockSize - 用于角点检测的邻域大小，就是上面提到的窗口的尺寸
- ksize - 用于计算梯度图的Sobel算子的尺寸
- k - 用于计算角点响应函数的参数k，取值范围常在0.04~0.06之间

```
import cv2 as cv
from matplotlib import pyplot as plt
import numpy as np

# detector parameters
block_size = 3
sobel_size = 3
k = 0.06

image = cv.imread('Scenery.jpg')

print(image.shape)
height = image.shape[0]
width = image.shape[1]
channels = image.shape[2]
print("width: %s  height: %s  channels: %s"%(width, height, channels))
   
gray_img = cv.cvtColor(image, cv2.COLOR_BGR2GRAY)


# modify the data type setting to 32-bit floating point 
gray_img = np.float32(gray_img)

# detect the corners with appropriate values as input parameters
corners_img = cv.cornerHarris(gray_img, block_size, sobel_size, k)

# result is dilated for marking the corners, not necessary
kernel = cv2.getStructuringElement(cv2.MORPH_RECT,(3, 3))
dst = cv.dilate(corners_img, kernel)

# Threshold for an optimal value, marking the corners in Green
#image[corners_img>0.01*corners_img.max()] = [0,0,255]

for r in range(height):
        for c in range(width):
            pix=dst[r,c]
            if pix>0.05*dst.max():
               cv2.circle(image,(c,r),5,(0,0,255),0)

image = cv.cvtColor(image, cv2.COLOR_BGR2RGB)
plt.imshow(image)
plt.show()
```





