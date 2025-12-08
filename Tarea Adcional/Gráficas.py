#!/usr/bin/env python
# coding: utf-8

# In[192]:


from sympy import *
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt


# In[193]:


# ======================
#    CURVA Ω_H vs Qp
# ======================
# Declarar variables simbólicas:
OmegaH, Qp = sp.symbols('Omega_{H}, Q_{p}')

# Plantear ecuación implícita de Ω_H:
equ_1 = sp.Eq(1 / 2, 1 / (1 + OmegaH**2 * (1 / (Qp)**2 - 2) + OmegaH**4))

# Despejar Ω_H:
sol_OmegaH = sp.solve(equ_1, OmegaH)


# In[194]:


# Definir variable de tipo 'Función'
OmegaH = sp.Function ('Omega_{H}') 

# Asignar solución
OmegaH = sol_OmegaH[3]


# In[195]:


# Eje X 
x = np.linspace(0.01, 1, 1000)

# Eje Y
y = 1.18920711500272 * np.sqrt(0.707106781186547 + ((x**4 - 0.5 * x**2 + 0.125)**0.5 / x**2) - 0.353553390593274 / x**2)


# In[196]:


# Configuración para Graficar Curva
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y, 'b-')
ax.grid(visible=True, linestyle='--', alpha=0.7)

ax.set_xlabel('Qp')
ax.set_ylabel('Ω_H')
ax.set_title('Curva: Ω_H (Qp)')

plt.show


# In[197]:


# ======================
#    CURVA Ω_G vs Qp
# ======================
# Declarar variables simbólicas:
OmegaG, Qp = sp.symbols('Omega_{G}, Q_{p}')

# Plantear ecuación implícita de Ω_G:
equ_2 = sp.Eq(1, OmegaG**4 + (OmegaG / Qp)**2)

# Despejar Ω_G:
sol_OmegaG = sp.solve(equ_2, OmegaG)


# In[198]:


# Definir variable de tipo 'Función'
OmegaG = sp.Function ('Omega_{G}') 

# Asignar solución
OmegaG = sol_OmegaG[3]


# In[199]:


# Eje X 
x = np.linspace(0.01, 1, 1000)

# Eje Y
y = np.sqrt(2) * np.sqrt(np.sqrt(4 * x**4 + 1) / x**2 - 1 / x**2) / 2


# In[200]:


# Configuración para Graficar Curva
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y, 'r-')
ax.grid(visible=True, linestyle='--', alpha=0.7)

ax.set_xlabel('Qp')
ax.set_ylabel('Ω_G')
ax.set_title('Curva: Ω_G (Qp)')

plt.show


# In[201]:


# ======================
#    CURVA Mp vs Qp
# ======================
# Declarar variables simbólicas:
OmegaG, Mp, Qp = sp.symbols('Omega_{G}, M_{p}, Q_{p}')

# Despejar Ω_G de equ_2
sol_OmegaG = sp.solve(equ_2, OmegaG)

# Plantear ecuación de Mp:
equ_3 = sp.Eq(Mp, sp.rad(90) - sp.atan(sol_OmegaG[3] * Qp))


# In[202]:


# Definir variable de tipo 'Función'
Mp = sp.Function ('Omega_{G}') 

# Asignar equ_3
Mp = equ_3


# In[203]:


# Eje X 
x = np.linspace(0.01, 1, 1000)

# Eje Y
y = np.degrees(-np.atan(np.sqrt(2) * x * np.sqrt((np.sqrt(4 * x**4 + 1) - 1) / x**2) / 2) + np.pi / 2)


# In[204]:


# Configuración para Graficar Curva
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y, 'g-')
ax.grid(visible=True, linestyle='--', alpha=0.7)

ax.set_xlabel('Qp')
ax.set_ylabel('Mp')
ax.set_title('Curva: Mp (Qp)')

plt.show

