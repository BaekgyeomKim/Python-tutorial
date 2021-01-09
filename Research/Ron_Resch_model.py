import numpy as np
import matplotlib.pyplot as plt

## input variable
alpha_1 = np.arange(0*np.pi, np.pi,  np.pi/1800)
alpha_1 = np.append(alpha_1,np.pi)
alpha_2 = np.arange(2/3*np.pi , np.pi, np.pi/5400)
print(alpha_1.shape)
print(alpha_2.shape)

R = 25.
print(type(R))

y = np.sin(alpha_1)
theta = 2*np.arcsin(np.sqrt((1/8)*(1-np.cos(alpha_1))))
print(theta.shape)

eta = np.arccos((np.sqrt(2*(1-np.cos(alpha_2)))-np.sin(theta/2))/(np.sqrt(3)*np.cos(theta/2)))
print(eta.shape)

## coordinate of the point 24
x_24 = np.array([])
y_24 = np.array([])

if np.sqrt(3)*np.tan(theta/2) > np.cos(eta):
    x_24 = np.sqrt(3)*R*np.cos(eta)*np.cos(theta/2)
else:
    x_24 = R/2*(np.sqrt(3)*np.cos(eta)*np.cos(theta/2)+3*np.sin(theta/2))

print(x_24)

# if np.sqrt(3)*np.tan(theta/2) > np.cos(eta):
#     y_24 = R*np.cos(eta)*np.cos(theta/2)
# else:
#     y_24 = R/2*(np.cos(eta)*np.cos(theta/2)+np.sqrt(3)*np.sin(theta/2))

z_24 = R*np.cos(theta/2)*np.sin(eta)


##


## plot
plt.plot(alpha_1, theta)
plt.show()
