import numpy as np
import math
from math import cos
from math import sin
from math import sqrt
from math import pi
import matplotlib.pyplot as plt

#condition implicite pour la suite: x<=L, et y<=l
c = 340*100   #Célérité de la lumière en cm par sec
t = 1  #durée d'étude pour les signaux sinusoïdaux
fe = 5000 #fréquence d'échantillonage

def echantillonnagecos(t1, f): #création de l'échelle des temps qui servira aux entrées sinusoïdales en fonction de la période d'analyse voulue. il faut que t1 coincide avec le t au dessus et que f coincide avec fe au dessus.
    T1 = np.linspace(0, t1, int(t1*f))         #échelle des temps qui servira pour les entrées sinusoïdales
    return T1
T = echantillonnagecos(1, fe)

# Création des microphones virtuelles, disposés selon le périmètre d'un cercle de rayon r à intervalles de corde régulier
def distance(L, l, r, x, y, F):
    X = np.zeros(32)
    Y = np.zeros(32)
    R = np.zeros(32)
    for i in range(32):
        X[i] = L/2 + r * cos(i * 2 * pi/32)
        Y[i] = l/2 + r * sin(i * 2 * pi/32)
        R[i] = sqrt(((x - X[i])**2)+((y - Y[i])**2) + F**2)
    return R

L = 200 # Length of the anechoïc chamber (cm)
l = 200 # Width of the anechoïc chamber (cm)
r = 150 # Radius of the sound rotation (cm)
x = 90 # Absciss of the measurement point on the focal plane (cm)
y = 120 # Ordinate of the measuremant point on the focal plane (cm)
F = 200 # Height of the focal plane (cm)
A = 5 # Amplitude of the sound signal (cm)
w = 6280 # Signal pulse

 # Création d'une liste de fonctions regroupant les signaux envoyées par les micros ou les entrées
def createemissioncos(L, l, r, x, y, F, A, w):
    E = np.zeros((32, 5000))  # liste de liste de fonctions
    for i in range(32):
        for j in range(5000):
            E[i][j] = (A * cos(w * T[j]))
    return E

W = createemissioncos(L, l, r, x, y, F, A, w)

# Cas de l'émission simultanée (sans rotation donc sans effet doppler)

# Calcul du signal reçu par un point M (x, y, F) du plan focal, comme superposition des signaux reçus
def signalcos(L, l, F, r, x, y, A, w, tm, W):     #tm est le temps de mesure compris entre 0 et 1
    S = np.zeros(32)
    D = distance(L, l, r, x, y, F)
    for i in range(32):
        if ((tm - (D[i] / c)) * 5000) >= 0 and ((tm - (D[i] / c)) * 5000) < 5000:
            S[i] = W[i][int((tm - D[i] / c) * 5000)]
        else:
            S[i] = (0)  #le signal n'est toujours pas arrivé au point d'intérêt
    return sum(S)

# Calcul du signal reçu par la droite d'équation x=110 du plan focal à l'instant tm
def cartographieligne(L, l, r, F, A, w, tm, W):  #on choisit de faire l'image de la ligne au milieu du plan focal en parcourant les ordonnées à l'instant tm
    K = np.zeros(200)
    for j in range(200):
        K[j] = signalcos(L, l, F, r, 110, j, A, w, tm, W)
    return K

# Calcul du signal reçu par la droite d'équation x=110 du plan focal en fonction du temps
def cartographielignetemporelle(L, l, r, F, A, w, W):
    G = np.zeros((5000,200))
    for i in range(len(T)):
        G[i] = cartographieligne(L, l, r, F, A, w, T[i], W)
    return G

# Calcul du signal reçu sur l'intégralité du plan focal à l'instant tm
def cartographie2Dinstant(L, l, r, F, A, w, tm, W):
    K = np.zeros((200, 200))
    for i in range(200):
        for j in range(200):
            K[i][j] = signalcos(L, l, F, r, i, j, A, w, tm, W)
    return K

# Calcul du signal reçu sur l'intégralité du plan focal à l'instant tm
def cartographie2Dtemporelle(L , l, r, F, A, w, W):
    N = np.zeros((200, 200, 5000))
    for k in range(5000):
        for i in range(200):
            for j in range(200):
                N[i][j][k] = signalcos(L, l, F, r, i, j, A, w, T[k], W)
    return N
            
# Cas de l'émission avec rotation en activant les micros les uns après les autres

#création du signal sinusoïdal tournant
def createemissionrota(L, l, r, x, y, F, A, w):
    G=createemissioncos(L, l, r, x, y, F, A, w)
    A = np.zeros((32, 5000))
    G[0][156 : 5000] = A[0][156 : 5000]
    G[31][0 : 4836] = A[31][0 : 4836]
    for i in range(1, 31):
        G[i][0 : i * 156] = A[i][0 : i * 156]
        G[i][(i + 1) * 156 : 5000] = A[i][(i + 1) * 156 : 5000]
    return G
    


def createemissionrota2(L, l, r, x, y, F, A, w): #pour  l'émission à 5 tours par seconde
    E1 = np.zeros((32, 1000))  # liste de liste de fonctions, regroupe les signaux envoyées par les micros ou les entrées
    for i in range(32):
        for j in range(1000):
            E1[i][j] = (A * cos(w * T[j]))
    A = np.zeros((32, 1000))
    E1[0][31 : 1000] = A[0][31 : 1000]
    E1[31][0 : 969] = A[31][0 : 969]
    for i in range(1, 31):
        E1[i][0 : i * 31]=A[i][0 : i * 31]
        E1[i][(i + 1) * 31 : 1000]=A[i][(i + 1) * 31 : 1000]
    return np.concatenate((E1, E1, E1, E1, E1), axis=1)

V=createemissionrota2(L, l, r, x, y, F, A, w)
    
    
# Calcul du signal reçu par un point M (x, y, F) du plan focal, comme superposition des signaux reçus, dans le cas d'un rotation
def signalcosrota(L, l, F, r, x, y, A, w, tm, V):     #tm est le temps de mesure compris entre 0 et 1
    S = np.zeros(32)
    D = distance(L, l, r, x, y, F)
    for i in range(32):
        if ((tm - (D[i] / c)) * 5000) >= 0 and ((tm - (D[i] / c)) * 5000) < 5000: #Si le signal est arrivé sur la plan focal
            S[i] = V[i][int((tm - D[i] / c) * 5000)]
        else:
            S[i] = (0) 
    return sum(S)


# Calcul du signal reçu par la droite d'équation x=110 du plan focal à l'instant tm, dans le cas d'un rotation
def cartographielignerota(L, l, r, F, A, w, tm, V):  #on choisit de faire l'image de la ligne au milieu du plan focal en parcourant les ordonnées à l'instant tm
    K = np.zeros(200)
    for j in range(200):
        K[j]=signalcosrota(L, l, F, r, 110, j, A, w, tm, V)
    return K

# Calcul du signal reçu par la droite d'équation x=110 du plan focal en fonction du temps, dans le cas d'un rotation
def cartographielignetemporellerota(L, l, r, F, A, w, V):
    G = np.zeros((5000, 200))
    for i in range(len(T)):
        G[i] = cartographielignerota(L, l, r, F, A, w, T[i], V)
    return G

# Calcul du signal reçu suur l'intégralité du plan focal à l'instant tm dans le cas d'un rotation
def cartographie2Dinstantrota(L, l, r, F, A, w, tm, V):
    K = np.zeros((200,200))
    for i in range(200):
        for j in range(200):
            K[i][j] = signalcosrota(L, l, F, r, i, j, A, w, tm, V)
    return K
# Calcul du signal reçu suur l'intégralité du plan focal à l'instant tm dans le cas d'un rotation

def cartographie2Dtemporellerota(L, l, r, F, A, w, V):
    N = np.zeros((5000, 200, 200))
    for k in range(5000):
        for i in range(200):
            for j in range(200):
                N[k][i][j] = signalcosrota(L, l, F, r, i, j, A, w, T[k], V)
    return N

# Retour du signal sur le plan où se trouvent les micros à l'instant tm
def retourmicroinstant(L, l, F, r, x, y, A, w, tm, V): 
    D = distance(L, l, r, x, y, F)
    K = np.zeros(32)
    for i in range(32):
        if ((tm - 2 * (D[i] / c)) * 5000) >= 0 and ((tm - 2 * (D[i] / c)) * 5000) < 5000:
            K[i] = signalcosrota(L, l, F, r, x, y, A, w, tm - D[i] / c, V)
        else:
            K[i] = 0
    return K
    
# Retour du signal sur le plan où se trouvent les micros en fonction du temops
def retourmicrotemporel(L, l, F, r, x, y, A, w, V):
    D = distance(L, l, r, x, y, F)
    K = np.zeros(5000, 32)
    for i in range(5000):
        for j in range(32):
            if ((T[i] -2 *(D[j] / c)) * 5000) >=0 and ((T[i] - 2 * (D[j] / c)) * 5000) < 5000:
                K[i][j] = signalcosrota(L, l, F, r, x, y, A, w, T[i] - D[j] / c, V)
            else:
                K[i][j] = 0
    return K

# Consutrction du profil final des 32 micros
def constructionsignal(L, l, F, r, x, y, A, w, V):
    L = np.zeros(5000, 1)
    L[0:155] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[0:155][0]
    L[156:310] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[156:310][1]
    L[311:465] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[311:465][2]
    L[466:620] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[466:620][3]
    L[621:775] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[621:775][4]
    L[776:930] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[776:930][5]
    L[931:1085] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[931:1085][6]
    L[1086:1240] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1086:1240][7]
    L[1241:1395] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1241:1395][8]
    L[1396:1550] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1396:1550][9]
    L[1551:1705] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1551:1705][10]
    L[1706:1860] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1706:1860][11]
    L[1861:2015] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[1861:2015][12]
    L[2016:2170] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2016:2170][13]
    L[2171:2325] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2171:2325][14]
    L[2326:2480] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2326:2480][15]
    L[2481:2635] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2481:2635][16]
    L[2636:2791] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2636:2791][17]
    L[2792:2946] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2792:2946][18]
    L[2947:3101] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[2947:3101][19]
    L[3102:3256] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3102:3256][20]
    L[3257:3411] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3257:3411][21]
    L[3412:3566] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3412:3566][22]
    L[3567:3721] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3567:3721][23]
    L[3722:3876] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3722:3876][24]
    L[3877:4031] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[3877:4031][25]
    L[4032:4186] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4032:4186][26]
    L[4187:4341] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4187:4341][27]
    L[4342:4496] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4342:4496][28]
    L[4497:4651] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4497:4651][29]
    L[4652:4806] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4652:4806][30]
    L[4807:4999] = retourmicrotemporel(L, l, F, r, x, y, A, w, V)[4807:4999][31]
    return L

# Cartographie finale contenant les rotations et les retours de son
def cartographiefinale(L, l, r, F, A, w, V):
    H = np.zeros((5000, 200, 200))
    for i in range(5000):
        for j in range(200):
            for k in range(200):
                H[:][j][k]=constructionsignal(L, l, F, r, j, k, A, w, V)
    return H
