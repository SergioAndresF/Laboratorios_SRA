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
# =====================
#   ESPECIFICACIONES
# =====================
Avf_dB    = 20       #  Ganancia de Lazo Cerrado en dB
Avf_veces = 10       #  Ganancia de Lazo Cerrado en veces
Mp        = 65       #  Margen de Fase para MPM
Qp        = 0.707    #  Factor de Calidad del Polo para MPM
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

# Para lograr Maxima Planicidad de Modulo (Mp = 65)
# Se debe cumplir que: fg / fpx = 0.4663 
# siendo: fg  -> Punto Critico
#         fpx -> Polo del AO2
fg_fpx = np.tan(np.radians(180 - 90 - Mp))

# Entonces por GBW:
fg  = np.sqrt(fT * Ad0_veces * fg_fpx)
fpx = fg / fg_fpx 

# Para ubicar fpx en el punto calculado
# Se debe cumplir que R2 = 1.159 * R1:
R2_R1 = (fT / fpx) - 1

# Entonces:
R1   = 1E3            #  R1 de la Red de Realimentacion de AO2
R2   = R2_R1 * R1     #  R2 de la Red de Realimentacion de AO2
Avf2 = R2_R1 + 1      #  Ganancia de Lazo Cerrado Ideal de AO2

# Calculos Sistema Compuesto
To   = K * Ad0_veces * Avf2      #  Ganancia de Lazo
fp   = np.sqrt(To * fp1 * fpx)   #  Frecuencia de Polo 
Dp   = fpx / fp1                 #  Distancia entre Polos 
Qp_s = np.sqrt(To / Dp)          #  Factor de Calidad del Polo 
fH   = 10**(3 / 20) * fg         #  Ancho de Banda a -3dB
# ===================================
#   RESPUESTA AL ESCALON UNITARIO
# ===================================
# De acuerdo a lo obtenido en MATLAB:
RiseTime      = 1.0539e-06
TransientTime = 3.0258e-06
SettlingTime  = 3.0258e-06
SettlingMin   = 9.1140
SettlingMax   = 10.4883
Overshoot     = 4.8880
Undershoot    = 0
Peak          = 10.4883
PeakTime      = 2.1931e-06
# ===================================
#   ESTIMACION DEL MARGEN DE FASE
# ===================================
# Calculamos Zeta (Coeficiente de Amortiguamiento):
Zeta = - np.log(Overshoot / 100) / np.sqrt(((np.pi)**2 + (np.log(Overshoot / 100))**2))

# Calculamos Mp:
Mp_step = np.degrees(np.arctan(2 * Zeta / (np.sqrt(np.sqrt(1 + 4 * Zeta**4) - 2 * Zeta**2))))