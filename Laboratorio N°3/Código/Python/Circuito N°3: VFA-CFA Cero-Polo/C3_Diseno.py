# =============
#   PAQUETES
# =============
import numpy as np        
import sympy as sp
from scipy.signal import ZerosPolesGain, TransferFunction, bode
import pandas as pd
import matplotlib.pyplot as plt
# ==========================
#   CARACTERISTICAS LM324
# ==========================
Ad0_dB    = 100       #  Ganancia Diferencial en dB
Ad0_veces = 100E3     #  Ganancia Diferencial en veces
fT        = 1e6       #  Frecuencia de Cruce
fp1       = 10        #  1er Polo
fp2       = 5.06E6    #  2do Polo
# ===========================
#   CARACTERISTICAS LM6181
# ===========================
RT   = 2.37E6      #  Transresistencia
CT   = 4.8E-12     #  Transcapacitancia
fpx  = 33.34E6     #  Polo del CFA sin Compensador
Avf2 = 20          #  Ganancia del CFA sin Compensador
R2   = 994.6       #  R2 de AO2 sin Compensador
R1   = 52.3        #  R1 de AO2 sin Compensador
# =====================
#   ESPECIFICACIONES
# =====================
Avf_dB    = 20       #  Ganancia de Lazo Cerrado en dB
Avf_veces = 10       #  Ganancia de Lazo Cerrado en veces
# ============================
#   COMPENSADOR CERO - POLO
# ============================
fzc  = 5.06E6        #  Cero, cancela 2do Polo del VFA
fpc  = 2 * 5.06E6    #  Polo, a una Octava del Cero
Comp = fzc / fpc     #  Atenuacion producida por el Compensador

# Para ubicar adecuadamente fzc y fpc
# Se debe cumplir que: R3 = R4
# siendo: R3 -> Resistencia en Paralelo al Capacitor
#         R4 -> Resistencia a Masa
R3 = 1E3
R4 = 1E3
# Y para hallar C:
C = 1 / (2 * np.pi * fzc * R3)    # Capacitor 

# ===========================
#   AMPLIFICADOR COMPUESTO
# ===========================
# =================
#   DESARROLLOS
# =================
# Para lograr una Ganancia de Lazo Cerrado de 20 [dB]
# Se debe cumplir que: Rf = 9 * Ri
Ri = 1E3              #  Resistencia Ri
Rf = 9 * Ri           #  Resistencia Rf
# Entonces:
K  = Ri / (Ri + Rf)   #  Cantidad de Realimentacion

# Para compensar la Atenuacion
Avf2 = Avf2 / Comp
# Entonces:
R2 = 994.6            #  Manteniendo ubicacion de fpx
R1 = R2 / (Avf2 - 1)  #  Nuevo valor de R1 de AO2

# Hallar fg por GBW:
fg = Ad0_veces * Avf2 * Comp * fp1 * K
# El Mp queda:
Mp_s = (np.radians(180 - 90) - np.arctan(fg / fpx) - np.arctan(fg / fpc)) * 180 / np.pi

# Calculos Sistema Compuesto
To   = K * Ad0_veces * Avf2 * Comp   #  Ganancia de Lazo
fp   = np.sqrt(To * fp1 * fpc)       #  Frecuencia de Polo
Dp   = fpc / fp1                     #  Distancia entre Polos
Qp_s = np.sqrt(To / Dp)              #  Factor de Calidad del Polo
fH   = 10**(3 / 20) * fg             #  Ancho de Banda a -3dB
# ===================================
#   RESPUESTA AL ESCALON UNITARIO
# ===================================
# De acuerdo a lo obtenido en LTspice:
Peak          = 10.063728
Overshoot     = (Peak - 10) / 10 * 100
# ===================================
#   RESPUESTA AL ESCALON UNITARIO
# ===================================
# Calculamos Zeta (Coeficiente de Amortiguamiento):
Zeta = - np.log(Overshoot / 100) / np.sqrt(((np.pi)**2 + (np.log(Overshoot / 100))**2))

# Calculamos Mp:
Mp_step = np.degrees(np.arctan(2 * Zeta / (np.sqrt(np.sqrt(1 + 4 * Zeta**4) - 2 * Zeta**2))))