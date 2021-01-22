import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pandas as pd
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

## colormap setting
viridis = cm.get_cmap('viridis', 18)
viridis6 = cm.get_cmap('viridis', 6)
plasma = cm.get_cmap('plasma', 6)
cividis = cm.get_cmap('cividis', 6)

## input variables
alpha_1 = np.arange(0*np.pi, np.pi,  np.pi/1800)
alpha_1 = np.append(alpha_1,np.pi)
alpha_1d = np.rad2deg(alpha_1)
alpha_2 = np.arange(2/3*np.pi , np.pi, np.pi/5400)
alpha_2d = np.rad2deg(alpha_2)
# print(alpha_1.shape)
# print(alpha_2.shape)

alpha_1d_data = pd.DataFrame(alpha_1d)
# alpha_1d_data.to_csv("alpha_1d.csv", header=False, index=False, encoding='utf-8')

R = 25.
# print(type(R))

y = np.sin(alpha_1)
theta = 2*np.arcsin(np.sqrt((1/8)*(1-np.cos(alpha_1))))
theta_d = np.rad2deg(theta)
# print(theta.shape)

eta = np.arccos((np.sqrt(2*(1-np.cos(alpha_2)))-np.sin(theta/2))/(np.sqrt(3)*np.cos(theta/2)))
eta_d = np.rad2deg(eta)
# print(eta.shape)

## coordinates of the point 24
X_24 = np.array([])
Y_24 = np.array([])

for i in range(1801):
    if np.sqrt(3)*np.tan(theta[i]/2) > np.cos(eta[i]):
        x_24 = np.sqrt(3)*R*np.cos(eta[i])*np.cos(theta[i]/2)
        X_24 = np.append(X_24, x_24)
    else:
        x_24 = R/2*(np.sqrt(3)*np.cos(eta[i])*np.cos(theta[i]/2)+3*np.sin(theta[i]/2))
        X_24 = np.append(X_24, x_24)

# print(X_24)
# print(X_24.shape)

for i in range(1801):
    if np.sqrt(3)*np.tan(theta[i]/2) > np.cos(eta[i]):
        y_24 = R*np.cos(eta[i])*np.cos(theta[i]/2)
        Y_24 = np.append(Y_24, y_24)
    else:
        y_24 = R/2*(np.cos(eta[i])*np.cos(theta[i]/2)+np.sqrt(3)*np.sin(theta[i]/2))
        Y_24 = np.append(Y_24, y_24)
# print(Y_24)
# print(Y_24.shape)

Z_24 = R*np.cos(theta/2)*np.sin(eta)

P_24 = np.transpose(np.array([X_24, Y_24, Z_24]))
# print(P_24)

## coordinates of the point 14

X_14 = R/2*(np.sqrt(3)*np.cos(eta)*np.cos(theta/2)+np.sin(theta/2))
Y_14 = R/2*np.abs(np.sqrt(3)*np.sin(theta/2)-np.cos(eta)*np.cos(theta/2))
Z_14 = Z_24

P_14 = np.transpose(np.array([X_14, Y_14, Z_14]))

## coordinates of the point 22

X_22 = np.array([])
Y_22 = np.array([])

for i in range(1801):
    if np.sqrt(3)*np.tan(theta[i]/2) > np.cos(eta[i]):
        x_22 = R*np.sin(theta[i]/2)
        X_22 = np.append(X_22, x_22)
    else:
        x_22 = R/2*(np.sqrt(3)*np.cos(eta[i])*np.cos(theta[i]/2)-np.sin(theta[i]/2))
        X_22 = np.append(X_22, x_22)

# print("X_22 shape is", X_22.shape)

for i in range(1801):
    if np.sqrt(3)*np.tan(theta[i]/2) > np.cos(eta[i]):
        y_22 = R*np.cos(eta[i])*np.cos(theta[i]/2)
        Y_22 = np.append(Y_22, y_22)
    else:
        y_22 = R/2*(np.cos(eta[i])*np.cos(theta[i]/2)+np.sqrt(3)*np.sin(theta[i]/2))
        Y_22 = np.append(Y_22, y_22)

# print("Y_22 shape is", Y_22.shape)

Z_22 = Z_24

# print("Z_22 shape is", Z_22.shape)

P_22 = np.transpose(np.array([X_22, Y_22, Z_22]))

## coordinates of the point 13

X_13 = np.zeros([1801,])
# print("X_13 shape is", X_13.shape)

Y_13 = np.zeros([1801,])
Z_13 = np.zeros([1801,])

P_13 = np.transpose(np.array([X_13, Y_13, Z_13]))

## coordinates of the point 15

X_15 = R*(np.sqrt(3)*np.cos(eta)*np.cos(theta/2)+np.sin(theta/2))
Y_15 = np.zeros([1801,])
Z_15 = np.zeros([1801,])

P_15 = np.transpose(np.array([X_15, Y_15, Z_15]))

## coordinates of the point 32

X_32 = X_15/2
Y_32 = np.sqrt(3)/2*X_15
Z_32 = np.zeros([1801,])

P_32 = np.transpose(np.array([X_32, Y_32, Z_32]))

## coordinates of the point 23

X_23 = X_14
Y_23 = np.sqrt(3)/6*X_15
Z_23 = -np.sqrt(np.power(np.sqrt(3)/3*2*R, 2)-np.power(X_23, 2)-np.power(Y_23, 2))

P_23 = np.transpose(np.array([X_23, Y_23, Z_23]))

## distance from point 14 to point 22

d = np.sqrt(np.power(X_14 - X_22, 2) + np.power(Y_14 - Y_22, 2))

## np.array를 통해서 [x,y] 좌표가 묶인 array를 생성 후 plt.Polygon 함수를 통해 삼각형 그리기.

## estimating ellipsoid equation

A = np.sqrt(np.power(X_23, 2) + np.power(Y_23, 2))
B = A
C = np.power(Z_14, 2)*np.power(A, 2)/(X_15*X_23 + 2*Y_14*Y_23 - (np.power(X_15, 2))/4 - (np.power(Y_14, 2)))

u = np.linspace(0.0, 2.0*np.pi, 60)
v = np.linspace(0.0, 1.0/2.0*np.pi, 30)
t = np.linspace(0.0, np.pi, 30)

xx1 = A[316, ]*np.cos(t)
yy1 = C[316, ]*np.sin(t)
xx2 = A[332, ]*np.cos(t)
yy2 = C[332, ]*np.sin(t)
xx3 = A[341, ]*np.cos(t)
yy3 = C[341, ]*np.sin(t)
xx4 = A[369, ]*np.cos(t)
yy4 = C[369, ]*np.sin(t)
xx5 = A[453, ]*np.cos(t)
yy5 = C[453, ]*np.sin(t)
xx6 = A[683, ]*np.cos(t)
yy6 = C[683, ]*np.sin(t)

x_zero = np.zeros((30,1))

## radius of curvature at point 13 and 24

def curvature(a,c,x,z):
    kappa = 1/(np.power(a, 2)*np.power(c, 2))*np.power(np.power(x, 2)/np.power(a, 4) + np.power(z, 2)/np.power(c, 4), -3./2.)
    return kappa

k_p24 = np.array([])
k_p13 = np.array([])

for i in range(1801):
    k24 = curvature(A[i, ], C[i, ], X_24[i, ], Z_24[i, ])
    k_p24 = np.append(k_p24, k24)

radius_p24 = 1./k_p24

# print(k_p24.shape)
# print(k_p13.shape)
# print(X_13.shape)
# print(alpha_1d.shape)

## experimental result
## deg 31.6, 33.2, 34.1, 36.9, 45.3, 68.3
# print(alpha_1d[316, ])
# print(alpha_1d[332, ])
# print(alpha_1d[341, ])
# print(alpha_1d[369, ])
# print(alpha_1d[453, ])
# print(alpha_1d[683, ])

x_24e = np.array([27.15, 27.97, 28.12, 28.95, 30.37, 33.68])
y_24e = np.array([18.21, 18.56, 18.65, 18.99, 19.74, 21.25])
z_24e = np.array([2.27, 2.52, 2.7, 3.7, 5.73, 7.64])
# print(x_24e)
# print(z_24e)

## plot
# plt.figure(1)
# plt.plot(alpha_1d, X_24)
# plt.plot(alpha_1d, X_23)
# plt.xlim(0,180)
# plt.ylim(20, 40)
# plt.grid(True)
# plt.xlabel(r'Dihedral angle ($\alpha_1$) [$\degree$]')
# plt.ylabel('X position [mm]')
# plt.legend(['Point 24', 'Point 23'])
#
# plt.figure(2)
# plt.plot(alpha_1d, Z_24,c='r')
# # plt.plot(alpha_1d, Z_23)
# plt.xlim(0,180)
# plt.ylim(0, 6)
# plt.grid(color='lightgray', linestyle='--')
# plt.xlabel(r'Dihedral angle ($\alpha_1$) [$\degree$]')
# plt.ylabel('Z position [mm]')
# plt.legend(['Point 24', 'Point 23'])
#
# plt.figure(3)
# triangles = [[0, 1, 2]]
# triang = np.array([])
# Triang = np.array([])
# for i in range (18):
#     triang = mtri.Triangulation([P_22[100*i, 0], P_24[100*i, 0], P_14[100*i, 0]],
#                                    [P_22[100*i, 1], P_24[100*i, 1], P_14[100*i, 1]], triangles)
#     Triang = np.append(Triang, triang)
#
# for i in range (18):
#     plt.triplot(Triang[i,], c=viridis(i/18))
#
# plt.grid(color='lightgray', linestyle='--')
# plt.xlabel('x position')
# plt.ylabel('y position')
#
# # Point 24, 22, 14의 위치를 원점이 Point 23이 기준이 되도록 이동시킨다.

plt.figure(4)
x = A[316, ]*np.outer(np.cos(u), np.sin(v))
y = B[316, ]*np.outer(np.sin(u), np.sin(v))
z = C[316, ]*np.outer(np.ones_like(u), np.cos(v))
ax = Axes3D(plt.figure(4))
ax.plot_surface(x, y, z, cmap='plasma', alpha=0.4)
ax.plot(x_zero, xx1, yy1, color='k')
ax.set_zlim(0, 60)
ax.set_xlabel('x position')
ax.set_ylabel('y position')
ax.set_zlabel('z position')

plt.figure(5)
x = A[683, ]*np.outer(np.cos(u), np.sin(v))
y = B[683, ]*np.outer(np.sin(u), np.sin(v))
z = C[683, ]*np.outer(np.ones_like(u), np.cos(v))
ax = Axes3D(plt.figure(5))
ax.plot_surface(x, y, z, cmap='plasma', alpha=0.4)
ax.plot(x_zero, xx6, yy6, color='k')
ax.set_zlim(0, 60)
ax.set_xlabel('x position')
ax.set_ylabel('y position')
ax.set_zlabel('z position')

plt.figure(6)
plt.plot(xx1, yy1, c='k')
# plt.plot(xx2, yy2, c='k')
# plt.plot(xx3, yy3, c='k')
# plt.plot(xx4, yy4, c='k')
# plt.plot(xx5, yy5, c='k')
# plt.plot(xx6, yy6, c='k')
plt.xlabel("x position", fontsize=12)
plt.ylabel("z position", fontsize=12)
plt.xlim(-30, 30)
plt.ylim(0, 60)
plt.grid(color='lightgray', linestyle='--')

plt.figure(7)
ax = Axes3D(plt.figure(7))
ax.scatter(X_24[316, ], Y_24[316, ], Z_24[316, ], cmap='viridis', marker='o')
ax.scatter(X_24[332, ], Y_24[332, ], Z_24[332, ], cmap='viridis', marker='o')
ax.scatter(X_24[341, ], Y_24[341, ], Z_24[341, ], cmap='viridis', marker='o')
ax.scatter(X_24[369, ], Y_24[369, ], Z_24[369, ], cmap='viridis', marker='o')
ax.scatter(X_24[453, ], Y_24[453, ], Z_24[453, ], cmap='viridis', marker='o')
ax.scatter(X_24[683, ], Y_24[683, ], Z_24[683, ], cmap='viridis', marker='o')
for i in range(6):
    ax.scatter(x_24e[i, ], y_24e[i, ], z_24e[i, ], cmap='viridis', marker='^')
ax.set_zlim(0, 10)
ax.set_xlabel('x position')
ax.set_ylabel('y position')
ax.set_zlabel('z position')

plt.figure(8)
plt.plot(alpha_1d, k_p24, color='r')
plt.title(r'Curvature ($\kappa$) vs. Dihedral angle ($\alpha_1$)')
plt.xlabel(r'Dihedral angle, $\alpha_1$ [$\degree$]')
plt.ylabel(r'Curvature, $\kappa$ [1/mm]')
plt.grid(color='lightgray', linestyle='--')
plt.xlim(0, 180)
plt.ylim(0, 0.16)

plt.figure(9)
plt.plot(alpha_1d, radius_p24, color='r')
plt.title(r'Radius of curvature (R) vs. Dihedral angle ($\alpha_1$)')
plt.xlabel(r'Dihedral angle, $\alpha_1$ [$\degree$]')
plt.ylabel(r'Radius of curvature, R [mm]')
plt.grid(color='lightgray', linestyle='--')
plt.xlim(0, 180)
plt.ylim(0, 100)

plt.show()