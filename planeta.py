#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

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
        raiz=np.sqrt(x**2 + y**2)
        fx = (x / raiz) * (1- 2*self.alpha / raiz)
        fy = (y / raiz) * (1- 2*self.alpha / raiz)
        return [vx, vy, fx, fy]

    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
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
        pass

    def energia_total(self):
        '''
        Calcula la enérgía total del sistema en las condiciones actuales.
        '''
        pass
