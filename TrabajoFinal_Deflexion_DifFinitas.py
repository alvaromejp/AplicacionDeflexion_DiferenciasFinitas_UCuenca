# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 15:03:34 2021

@author: Alvaro Mejía

"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# Condiciones de Contorno

u_d = 0
u_c = 0

u_a = 0
u_b = 0

#    Caracterización del elemento estructural
# Resistencia a la compresión f'c en (MPa)
fprima_c = 50
# Modulo de elasticidad del agregado (GPa)
Ea = 110.5

Ec = 1.15*(Ea**(1/3))*(fprima_c**(1/2))

# Carga aplicada (Pa)
q = 33.6e3

# Constante de Poisson
sigma = 0.3

# Grosor (m)
delta_z = 0.1

# Modulo de elasticidad (kPa)

#E = Ec*1e6
E = 2e8
#Rigidez de flexión 

D = (E*delta_z**3)/(12*(1-sigma**2))

# Dimensiones de la losa

a = 2 # alto
b = 2 # bajo



# Tamaños de paso para x e y 
delta_xy = 0.5

max_iter = 100
error = 0.0001

res = (q*(delta_xy**2))/D
# Número de particiones de la malla
n = int(b/delta_xy)+1
m = int(a/delta_xy)+1

u = np.zeros((m,n))
for j in range(m): 
    u[j,0] = u_a
    u[j,-1] = u_b
for i in range(n): 
    u[0,i] = u_d
    u[-1,i] = u_c
    
p =(u_a + u_b + u_c + u_d)/4
for j in range(1,m):
    for i in range(1,n): 
        u[j,i] = p
i = 0
conv = 0

while i < max_iter and conv == 0: 
    i +=1
    t = u.copy()
    for j in range(1,m-1):
        for i in range(1,n-1): 
            u[j,i] = 0.25*(u[j,i-1] + u[j,i+1] + u[j+1,i]+u[j-1,i] -res)
    if np.linalg.norm(u-t, np.inf)/np.linalg.norm(u,np.inf) < error: 
        conv = 1
        
if conv ==1: 
    print(u)
    
    
print()
"""
Valores de Z
"""

z = np.zeros((m,n))
for j in range(m): 
    z[j,0] = u_a
    z[j,-1] = u_b
for i in range(n): 
    z[0,i] = u_d
    z[-1,i] = u_c
    
p =(u_a + u_b + u_c + u_d)/4
for j in range(1,m):
    for i in range(1,n): 
        z[j,i] = p
i = 0
conv = 0

while i<max_iter and conv == 0: 
    i +=1
    t = z.copy()
    for j in range(1,m-1):
        for i in range(1,n-1): 
            z[j,i] = 0.25*(z[j,i-1] + z[j,i+1] + z[j+1,i]+z[j-1,i] -(u[j,i])*(delta_xy**2))     
    if np.linalg.norm(z-t, np.inf)/np.linalg.norm(z,np.inf) < error: 
        conv = 1     
if conv ==1: 
    print(z)  
    
print()
print("Resistencia a la compresión f'c: ", fprima_c, "Mpa")
print("Módulo de elasticidad:",E*1e-6, "GPa")
print("Grosor de la losa:", delta_z, "m")

print()
print("Máxima deflexión (m):", round(abs(np.amax(z)),4))
print("Carga aplicada (N/m^2):", q)


fig = plt.figure(figsize=(8,6))
ax3d = plt.axes(projection="3d")

xdata = np.linspace(0,a,m)
ydata = np.linspace(0,b,n)
X,Y = np.meshgrid(xdata,ydata)
Z = z

ax3d = plt.axes(projection='3d')
ax3d.plot_surface(X, Y, Z,cmap='plasma')
#ax3d.set_title('Deflexión Ec ={0} GPa'.format(round(Ec)),fontsize=14,fontweight="bold")
ax3d.set_title('Deflexión en z',fontsize=14,fontweight="bold")
ax3d.set_xlabel('Longitud en metros')
ax3d.set_ylabel('Longitud en metros')
ax3d.set_zlabel("Deflexión en metros")
ax3d.set_zlim(0.10,0)


plt.show()

fig = plt.figure(figsize=(8,6))
plt.title("Representación gráfica de la\n deflexión utilizando diferencias finitas",fontsize=14,fontweight="bold")
plt.imshow(z,vmin=np.amin(z),vmax=np.amax(z),cmap='plasma')
plt.colorbar(shrink=.92)
plt.xticks(())
plt.yticks(())

plt.show()






























