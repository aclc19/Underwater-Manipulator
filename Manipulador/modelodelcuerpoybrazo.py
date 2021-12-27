# -*- coding: utf-8 -*-
"""TIF.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ODQRZseogHjwkU71axByI2-EvMisYNb7

# **MANIPULADOR 2GDL**
"""

import numpy as np
import math
import matplotlib.pyplot as plt

def d2r(grados):
  return grados * math.pi/180

l1 = 125 # Longitud del primer eslabon
l2 = 184 # Longitud del segundo eslabon

tita1_min = 64
tita1_max = 177
tita2_min = -125
tita2_max = 0

a1 = 110
a2 = 150

num = 20
theta11 = np.linspace(d2r(tita1_min),d2r(tita1_max),num)
theta22 = np.linspace(d2r(tita2_min), d2r(tita2_max),num)

fig = plt.figure ()
fig = plt.figure(figsize =((400)/60 ,(300+40)/60))
plt.ylim(300, -320)
plt.xlim( -300 ,300)
arg1 = theta11
arg2 = theta22 [0]+ theta11
px = (l1*np.cos(arg1) + l2*np.cos(arg2))
py = (l1*np.sin(arg1) + l2*np.sin(arg2))
x = px
y = py

arg1 = theta11[num -1]
arg2 = theta22+arg1
px = (l1*np.cos(arg1) + l2*np.cos(arg2))
py = (l1*np.sin(arg1) + l2*np.sin(arg2))
x = [*x,*px]
y = [*y,*py]

arg1 = theta11 [0]
arg2 = theta22+arg1
px = (l1*np.cos(arg1) + l2*np.cos(arg2))
py = (l1*np.sin(arg1) + l2*np.sin(arg2))
x = [*x,*px]
y = [*y,*py]

arg1 = theta11
arg2 = theta22[num -1]+ arg1
px = (l1*np.cos(arg1) + l2*np.cos(arg2))
py = (l1*np.sin(arg1) + l2*np.sin(arg2))
x = [*x,*px]
y = [*y,*py]

aux2 = x[2*num :4*num -1]
aux2 = aux2 [:: -1]
x_aux = [*x[0*num:2* num], *aux2]

aux2 = y[2*num :4*num -1]
aux2 = aux2 [:: -1]
y_aux = [*y[0*num:2* num], *aux2]

# EL espacio de trabajo
plt.fill(x_aux ,y_aux ,facecolor='lightgrey')

# Ploteo un brazo
JointNew = [100 ,30]
x1 = a1*np.cos(JointNew [0]*np.pi/180)
x2 = x1 + a2*np.cos(( JointNew [1]+ JointNew [0])*np.pi /180)
x = [0,x1,x2]
y1 = a1*np.sin(JointNew [0]*np.pi/180)
y2 = y1 + a2*np.sin(( JointNew [1]+ JointNew [0])*np.pi /180)
y = [0,y1,y2]

plt.plot(x,y, linewidth= 4)
#plt.plot(0,0,'o',markersize =10)
plt.plot(x1,y1,'o',markersize =10)
#plt.plot(x2,y2,'o',markersize =10)
#ax = plt.axes()

# El rov
theta = np.linspace(-np.pi, np.pi, 200)
plt.fill (150*np.sin(theta),150*np.cos(theta )-155,facecolor='blue', zorder = 10)
plt.xlabel("Eje x",fontsize =18)
plt.ylabel("Eje y",fontsize =18)
plt.grid(True , zorder = 5)
plt.text(-85, -130, 'Cuerpo ROV',fontsize =20, zorder = 10,color = 'white' )
plt.arrow( 0, 0, 120, 0, fc="k",ec="k", head_width =8,head_length =15, zorder = 10 )
plt.arrow( 0, 0, 0, 120, fc="k",ec="k", head_width =8,head_length =15, zorder = 10 )
plt.text (120, 20, 'x',fontsize =18)
plt.text(5, 140, 'y',fontsize =18)