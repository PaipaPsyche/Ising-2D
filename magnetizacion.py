# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 15:21:03 2018

@author: Nicolas
"""
import numpy as np
import matplotlib.pyplot as plt

#Recomendaciones:
#Dado que el codigo es tan ineficiente a nivel computacional
#se espera que se demore mucho corriendo
#una configuracion adecuada para correr sin esfuezo
#es dada por N<=22 ; t<=1500; n_it<20 ; n_Tintervalos<20


#tiempo esperado = 10 minutos

N=25 #cuadricula NxN
J=100000  #Parametro  J
t=10000  #numero de "ticks" de la simulación

kb=1 #Parámentro Kb
Tc=(2*J)/(kb*np.arcsinh(1))  #Temperatura critica
Tp=Tc*(0.4) #una teperatura mas baja que la temperatura critica
b=1/(kb*Tp) #Parametro beta

#Se define una matriz NxN y se llena aleatoriamente con 1 o -1
Mat=np.zeros(shape=(N,N))
def aleatorio():
    for i in range(N):
        for j in range(N):
            k=np.random.random()
            if(k>0.5):
                Mat[i,j]=1
            else:
                Mat[i,j]=-1
aleatorio()

#calculo de aporte de energia por pares de spin
def CalcE(Sa,Sb):
    return Sa*Sb

#calculo de la energia total del sistema con condiciones
#periodicas. Se suman las energias de todos los vinculos
#verticales y horizontales.

def CalcH(M):
    H=0
    #vECAL
    for j in range (N):
        for i in range (N):
           if(i==N-1):
               H=H+CalcE(M[i,j],M[0,j])
               
           else:
               H=H+CalcE(M[i,j],M[i+1,j])
              
    #horizontal
    for i in range (N):
        for j in range (N):
           if(j==N-1):
               H=H+CalcE(M[i,j],M[i,0])
              
           else:
               H=H+CalcE(M[i,j],M[i,j+1])
    return -J*H


#Se intenta cambiar un spin aleatorio del sistema
#se descarta si Delta>0 y p no cumple condiciones
#Se aprueba si Delta<=0 o si no, si p cumple condiciones
    
#Se retorna el valor total del cambio de energia
# es 0 si no se aprueba, y Delta si se aprueba.
#retorna la matriz nueva, cambiada o no.
def DeltaEAleatorio(Mat):
    MatRet=Mat.copy()
    MatPrueba=Mat.copy()
    #print(Mat)
    ret=0
    ir=int(np.random.random()*N)
    jr=int(np.random.random()*N)
    #print(ir,jr)
    Eini=CalcH(MatPrueba)
    #print(Mat)
    MatPrueba[ir,jr]=-MatPrueba[ir,jr]
    #print(Mat)
    Efin=CalcH(MatPrueba)
    Delta=Efin-Eini
    #print (Delta)
    if(Delta<=0):
       # print("a")
        ret=Delta
        MatRet=MatPrueba.copy()
        
    else:
        p=np.random.random()
        ee=np.exp(-b*Delta)
        if(p<ee):
            ret=Delta
            MatRet=MatPrueba.copy()
    return ret,MatRet

#Magnetizacion total de una matriz dada  
def magnetizacion(MAT):
    M=0
    for i in range (N):
        for j in range (N):
            M=M+MAT[i,j]
    
    return M/(N*N)

def EMfinal(t,MAT):
    Matp=MAT.copy()
    Mag=np.zeros([t,2])
    en=CalcH(Matp)
    for k in range (t):
        de,M=DeltaEAleatorio(Matp)
        Matp=M.copy()
        en=en+de
        magn=magnetizacion(Matp)
        Mag[k,0]=k
        Mag[k,1]=magn
	
    m=magnetizacion(Matp)
    return en,m,Mag

e,m,MAG=EMfinal(t,Mat)

plt.figure()
plt.plot(MAG[:,0],MAG[:,1])
plt.xlabel("Tiempo")
plt.ylabel("Magnetizacion")
plt.savefig("MagnetizacionIndividual.png")





