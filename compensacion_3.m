clc; clear; close all;

% ==========================
%   BLOQUE N°1: AO1 - VFA
% ==========================
Ad  = 100E3;              %  Ganancia Diferencial en veces
wp1 = 2 * pi * 10;        %  1er Polo
wp2 = 2 * pi * 5.06E6;    %  2do Polo

num_1 = Ad * wp1 * wp2;
den_1 = [1 (wp1 + wp2) (wp1 * wp2)];
G1    = tf(num_1, den_1);

% ============================
%   BLOQUE N°2: COMPENSADOR
% ============================
Comp = 0.5;                %  Atenuación 
wzc  = 2 * pi * 5.06E6;    %  Cero
wpc  = 2 * pi * 10.12E6;   %  Polo

num_2 = [(Comp * wpc) (Comp * wpc * wzc)];
den_2 = [wzc (wpc * wzc)];
G2    = tf(num_2, den_2);

% ==========================
%   BLOQUE N°3: AO2 - CFA
% ==========================
Av  = 40;                   %  Ganancia 
wpx = 2 * pi * 33.337E6;    %  Polo

num_3 = Av * wpx;
den_3 = [1 wpx];
G3    = tf(num_3, den_3);

% ===================
%   REALIMENTACIÓN
% ===================
K = 1/10;

% ===========================
%   AMPLIFICADOR COMPUESTO
% ===========================
GLA = G1 * G2 * G3;
GLC = feedback(GLA, K);

% ==================================
%   RESPUESTA AL ESCALÓN UNITARIO
% ==================================
step(GLC);
stepinfo(GLC);
