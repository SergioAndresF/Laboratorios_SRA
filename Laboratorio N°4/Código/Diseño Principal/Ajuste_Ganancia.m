clc; clear; close all;

%%
% ========================
%    AJUSTE DE GANANCIA 
% ========================
% Se tiene que:
Ganancia_dB    = 23;                        % Diferencia entre Curvas en dB
Ganancia_veces = 10^(Ganancia_dB / 20);     % Diferencia entre Curvas en veces
Atenuacion     = 1 / Ganancia_veces;        % Atenuacion buscada

% Entonces, si R1 = 100 [kOhm]:
R1 = 100E3;                   
R2 = Atenuacion * R1;         

% Resultados:
fprintf('      =================================\n');
fprintf('           Componentes Bloque N°3      \n');
fprintf('      =================================\n');

fprintf('Resistencias: R1 = %.0f [kOhm]; R2 = %.0f [kOhm]\n', R1 / 1E3, R2 / 1E3);
