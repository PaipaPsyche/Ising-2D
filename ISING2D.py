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

N=50 #cuadricula NxN
J=10  #Parametro  J
t=2000  #numero de "ticks" de la simulación
n_it=5  #numero iteraciones para sacae E y M promedio
n_Tintervalos=10    #intervalos de temperatura
kb=pow(10,-23)*1.38064852 #Parámentro Kb
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

#grafica la matriz y sus valores
def ploMat(MAT):
    plt.imshow(MAT)
#Magnetizacion total de una matriz dada  
def magnetizacion(MAT):
    M=0
    for i in range (N):
        for j in range (N):
            M=M+MAT[i,j]
    
    return M/(N*N)
#========================
E=[]
T=[]
e=CalcH(Mat)

E.append(e)
T.append(0)
MatS=Mat.copy()
#Se grafican el sistema inicial y el final de una iteracion
#Se grafica como varia la energia del sistema en el tiempo
FIGI=plt.figure(figsize=(10,10))
plt.title("Configuracion Inicial "+str(N)+"x"+str(N),fontsize=20)
ploMat(MatS)
FIGI.savefig("Cinicial.png")

for k in range (t):
    de,M=DeltaEAleatorio(MatS)
    MatS=M.copy()
    e=e+de
    E.append(e)
    T.append(k+1)
FIG=plt.figure(figsize=(10,10))
plt.plot(T,E)
plt.xlabel("Tiempo")
plt.ylabel("Energia")

bb=float("{0:.2f}".format(b))
plt.title("Ising con "+str(N*N)+" particulas -"+str(t)+" repeticiones (J="+str(J)+", B="+str(bb)+")",fontsize=20)
FIG.savefig("EvolucionEnergia.png")


FIGF=plt.figure(figsize=(10,10))
plt.title("Configuracion Final "+str(N)+"x"+str(N),fontsize=20)
ploMat(MatS)
FIGF.savefig("Cfinal.png")

#Calcula la Energia y la magnetizacion de el sistema final estable 
def EMfinal(t,MAT):
    Matp=MAT.copy()
    en=CalcH(Matp)
    for k in range (t):
        de,M=DeltaEAleatorio(Matp)
        Matp=M.copy()
        en=en+de
    m=magnetizacion(Matp)
    return e,m


#Calcula el promedio de magnetizacion y energia en los sitemas estables finales
def EMpromedio(t,n, MAT):
    Entot=0
    Mtot=0
    for h in range (n):
        e,m=EMfinal(t,MAT)
        Entot=Entot+e
        Mtot=Mtot+m
    return Entot/n,Mtot/n

       
#se barren ciertos intervalos de energia antes y despues de la
#temperatura critica. Se grafica magnetizacion Vs temperatura.Escala logaritmica.
expC=np.log10(Tc)
multiExp=np.linspace(int(0.1*expC),int(1.5*expC),n_Tintervalos)
multi=[]
for expe in multiExp:
    multi.append(pow(10,expe))
Mparr=[]
#se aplico escala logaritmica al rededor de la tmperatura critica
for temperature in multi:
    Tp=temperature
    Ep,Mp=EMpromedio(t,n_it,Mat)
    Mparr.append(Mp)
FIGM=plt.figure(figsize=(10,10))
plt.semilogx(basex=10)
plt.scatter(multi,Mparr)
plt.xlabel("Temperatura")
plt.ylabel("Magnetizacion")
plt.title("Magnetizacion Vs. Temperatura",fontsize=20)
plt.axvline(x=Tc,c='r')
FIGM.savefig("Magnetizacion.png")


  


