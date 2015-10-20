#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt



class Planeta(object):
    '''
    Complete el docstring.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual = condicion_inicial
        self.t_actual = 0.
        self.alpha = alpha

    def ecuacion_de_movimiento(self):
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.
        '''
        x, y, vx, vy = self.y_actual
        r = np.sqrt(x**2 + y**2)
        fx = -(1/r**2 + self.alpha/r**3)*(x/r)   # asumimos M,G,m = 1
        fy = -(1/r**2 + self.alpha/r**3)*(y/r)

        return [vx, vy, fx, fy]

    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
        vx, vy, ax, ay= self.ecuacion_de_movimiento()

        xn1 = self.y_actual[0] + dt * vx
        yn1 = self.y_actual[1] + dt * vy
        vxn1 = self.y_actual[2] + dt * ax
        vyn1 = self.y_actual[3] + dt * ay

        self.y_actual= [xn1, yn1, vxn1, vyn1]

        pass



    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        primero haremos runge kuta para tener X e Y, y luego para tener Vx y Vy
        '''
        vx1, vy1, ax1, ay1 = self.ecuacion_de_movimiento()
        K1 = vx1*dt, ax1*dt
        L1 = vy1*dt, ay1*dt
        self.y_actual[0] += K1[0]/2.0
        self.y_actual[1] += L1[0]/2.0

        vx2, vy2, ax2, ay2 = self.ecuacion_de_movimiento()
        K2 = vx2*dt, ax2*dt
        L2 = vy2*dt, ay2*dt
        self.y_actual[0] += K2[0]/2.0
        self.y_actual[1] += L2[0]/2.0

        vx3, vy3, ax3, ay3 = self.ecuacion_de_movimiento()
        K3 = vx3*dt, ax3*dt
        L3 = vy3*dt, ay3*dt
        self.y_actual[0] += K3[0]
        self.y_actual[1] += L3[0]

        vx4, vy4, ax4, ay4 = self.ecuacion_de_movimiento()
        K4 = vx4*dt, ax4*dt
        L4 = vy4*dt, ay4*dt

        '''
        ya hemos aumentado los parametros X e Y de las condiciones actuales, asi que ahora lo volveremos a operar para que calse
        con la formula de Runge-Kuta --> (1/6.0) * (K1[0] + 2*K2[0] + 2*K3[0] + K4[0]), analogo para L.
        Los elementos de velocidad de la condicion actual no los hemos modificado
        '''
        X_n1 = self.y_actual[0] - 2*K1[0]/6.0 - K2[0]/6.0 - 2*K3[0]/3.0 + K4[0]/6.0 #+ (1/6.0) * (K1[0] + 2*K2[0] + 2*K3[0] + K4[0])
        Y_n1 = self.y_actual[1] - 2*L1[0]/6.0 - L2[0]/6.0 - 2*L3[0]/3.0 + L4[0]/6.0 #+ (1/6.0) * (L1[0] + 2*L2[0] + 2*L3[0] + L4[0])

        ''' vx y vy cumplen con la forma original para Yn+1 de Runge-Kuta
        '''
        Vx_n1 = self.y_actual[2] + (1/6.0) * (K1[1] + 2*K2[1] + 2*K3[1] + K4[1])
        Vy_n1 = self.y_actual[3] + (1/6.0) * (L1[1] + 2*L2[1] + 2*L3[1] + L4[1])

        self.y_actual= [X_n1, Y_n1, Vx_n1, Vy_n1]


        pass



    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        [x, y, vx, vy] = self.y_actual
        vx1, vy1, ax1, ay1 = self.ecuacion_de_movimiento()

        x_n1 = x + dt*vx1 + ax1*(dt**2)/2.
        y_n1 = y + dt*vy1 + ay1*(dt**2)/2.

        self.y_actual = [x_n1, y_n1, vx1, vy1]

        vx2, vy2, ax2, ay2 = self.ecuacion_de_movimiento()

        vxf = vx + ax2 * dt/2.0 + ax1 * dt/2.0
        vyf = vy + ay2 * dt/2.0 + ay1 * dt/2.0

        self.y_actual = [x_n1, y_n1, vxf, vyf]

        pass

    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        Asumimos g=1, M=1 y m=1
        '''
        x, y, vx, vy = self.y_actual
        vx, vy, ax, ay = self.ecuacion_de_movimiento()

        r = np.sqrt(x**2+y**2)

        E = ((vx**2 + vy**2) * (0.5)) - 1/r + (self.alpha / (r**2))
        return E



        pass


'''
Ahora implementamos el codigo para plotear
'''

#planetaX = Planeta([10, 0, 0, 0.1], 10**(-2.871))
planetaX = Planeta([10, 0, 0, 0.4], 0)

largo = 50000
fin = 5000.0
t = np.linspace(0, fin, largo)

VectorX = []
VectorX = np.append(VectorX,planetaX.y_actual[0])
VectorY = []
VectorY = np.append(VectorY,planetaX.y_actual[1])
VectorE = []

dt = fin/largo

for i in range(0, largo):
    #planetaX.avanza_verlet(dt)
    #planetaX.avanza_euler(dt)
    planetaX.avanza_rk4(dt)
    VectorX = np.append(VectorX,planetaX.y_actual[0])
    VectorY = np.append(VectorY,planetaX.y_actual[1])
    VectorE = np.append(VectorE,planetaX.energia_total())
    planetaX.t_actual += dt


fig = plt.figure()

plt.plot(VectorX, VectorY, 'g')
plt.xlabel('x')
plt.ylabel('y')
plt.title("orbita para el metodo de RK4 con alfa= 0")
#plt.title("orbita para el metodo de RK4 con alfa = 10**(-2.871)")
plt.grid()

fig = plt.figure()

plt.plot(t, VectorE, 'g')
plt.xlabel('tiempo')
plt.ylabel('Energia')
plt.title("energia para el metodo de RK4")
plt.grid()

plt.show()
