clc; clear; close all;
%%
% ====================================
%   AJUSTE DE GANANCIA - BLOQUE N°1
% ====================================
R1_b1               = 3.14E3;
Peak_dB_obtenido_b1 = 14.396155;
Peak_dB_esperado_b1 = 2.2201632;
Atenuacion_b1       = 1 / 10^((Peak_dB_obtenido_b1 - Peak_dB_esperado_b1) / 20);

% Partiendo de que:
% ------------------------------------------
% EQ1: R1 = R1_1 // R1_2
% EQ2: Atenuacion = R1_2 / (R1_1 + R1_2)
% ------------------------------------------

% De EQ1 / EQ2:
R1_1_b1 = R1_b1 / Atenuacion_b1;

% Luego, de EQU1:
R1_2_b1 = - R1_b1 * R1_1_b1 / (R1_b1 - R1_1_b1);

% Resultados:
fprintf('            =================================\n');
fprintf('                 Divisor Resistivo Bloque N°1      \n');
fprintf('            =================================\n');

fprintf('Resistencias: R1_1 = %.2f [kΩ]; R1_2 = %.1f [kΩ]\n', R1_1_b1 / 1E3, R1_2_b1 / 1E3);

%%
% ====================================
%   AJUSTE DE GANANCIA - BLOQUE N°2
% ====================================
Peak_dB_obtenido_b2 = 14.609175;
Peak_dB_esperado_b2 = -258.99531E-3;
Atenuacion_b2       = 1 / 10^((Peak_dB_obtenido_b2 - Peak_dB_esperado_b2) / 20);

% Partiendo de que:
% ------------------------------------------
% EQ: Atenuacion = R1_2 / (R1_1 + R1_2)
% ------------------------------------------
R1_2_b2 = 1E3;

% De EQ1:
R1_1_b2 = (1 / Atenuacion_b2 - 1) * R1_2_b2;

% Resultados:
fprintf('            =================================\n');
fprintf('                 Divisor Resistivo Bloque N°2      \n');
fprintf('            =================================\n');

fprintf('Resistencias: R1_1 = %.2f [kΩ]; R1_2 = %.f [kΩ]\n', R1_1_b2 / 1E3, R1_2_b2 / 1E3);
