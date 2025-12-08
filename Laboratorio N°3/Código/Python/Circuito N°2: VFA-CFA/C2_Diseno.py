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
RT = 2.37E6      #  Transresistencia
CT = 4.8E-12     #  Transcapacitancia
# =====================
#   ESPECIFICACIONES
# =====================
Avf_dB    = 20       #  Ganancia de Lazo Cerrado en dB
Avf_veces = 10       #  Ganancia de Lazo Cerrado en veces
Mp        = 65       #  Margen de Fase para MPM
Qp        = 0.707    #  Factor de Calidad del Polo para MPM
fg        = 2E6      #  Ancho de Banda Potencial (Punto Critico)
# ===========================
#   AMPLIFICADOR COMPUESTO
# ===========================
# =================
#   DESARROLLOS
# =================
# Para lograr una Ganancia de Lazo Cerrado de 20 [dB]
# Se debe cumplir que: Rf = 9 * Ri
Ri = 1E3             #  Resistencia Ri
Rf = 9 * Ri          #  Resistencia Rf

# Entonces:
K  = Ri / (Ri + Rf)  #  Cantidad de Realimentacion

# Para lograr Maxima Planicidad de Modulo (Mp = 65)
# siendo: fpx -> Polo del AO2
fpx = fg / (np.tan(np.radians(180 - 90 - Mp) - np.arctan(fg / fp2)))

# Se sabe que: fpx = 1 / (2 * pi * CT * R2)
# Entonces:
R2 = 1 / (2 * np.pi * CT * fpx)     #  R2 de Red de AO2

# Y por GBW:
R1   = R2 / ((fg / Ad0_veces) - 1)  #  R1 de Red de AO2
Avf2 = 1 + R2 / R1                  #  Ganancia Lazo Cerrado AO2

# Calculos Sistema Compuesto
To   = K * Ad0_veces * Avf2       #  Ganancia de Lazo
fp   = np.sqrt(To * fp1 * fp2)    #  Frecuencia de Polo
Dp   = fp2 / fp1                  #  Distancia entre Polos
Qp_s = np.sqrt(To / Dp)           #  Factor de Calidad del Polo
fH   = 10**(3 / 20) * fg          #  Ancho de Banda a -3dB
# ===================================
#   RESPUESTA AL ESCALON UNITARIO
# ===================================
# De acuerdo a lo obtenido en MATLAB:
RiseTime      = 1.1462e-07
TransientTime = 2.9591e-07
SettlingTime  = 2.9591e-07
SettlingMin   = 9.0603
SettlingMax   = 10.3004
Overshoot     = 3.0041
Undershoot    = 0
Peak          = 10.3004
PeakTime      = 2.4275e-07
# ===================================
#   ESTIMACION DEL MARGEN DE FASE
# ===================================
# Calculamos Zeta (Coeficiente de Amortiguamiento):
Zeta = - np.log(Overshoot / 100) / np.sqrt(((np.pi)**2 + (np.log(Overshoot / 100))**2))
# Calculamos Mp:
Mp_step = np.degrees(np.arctan(2 * Zeta / (np.sqrt(np.sqrt(1 + 4 * Zeta**4) - 2 * Zeta**2))))