clc; clear; close all;
%%
% ======================
%    ESPECIFICACIONES
% ======================
fp1 = 800;        % Frecuencia inferior Pasabanda
fp2 = 1250;       % Frecuencia superior Pasabanda
fs1 = 200;        % Frecuencia inferior Rechazabanda
fs2 = 5000;       % Frecuencia superior Rechazabanda

Ap = 0.25;        % Atenuación en la Banda de Paso
As = 30;          % Atenuación en la Banda de Rechazo

%%
% ==================
%    DESARROLLOS
% ==================
% Convertir las frecuencias de [Hz] a [rad/seg]
Wp = 2 * pi * [fp1 fp2];
Ws = 2 * pi * [fs1 fs2];

% Calcular Orden Mínimo del Filtro (n) y Frecuencia Normalizada (Wc)
[n, Wc] = cheb1ord(Wp, Ws, Ap, As, 's');   
fprintf('Orden mínimo requerido: n = %d\n', n);

% Diseñar Filtro de Chebyshev Tipo 1
[num, den] = cheby1(n, Ap, Wc, 'bandpass', 's');  
Filtro     = tf(num, den);    
fprintf('La Función de Transferencia obtenida es:\n');

% Descomponer la Función de Transferencia en Bicuadráticas
[sos, g] = tf2sos(num,den);

PasaBajo = tf(2 * g * sos(1,1:3),sos(1,4:6));
PasaAlto = tf(1 / 2 * sos(2,1:3),sos(2,4:6));

%%
% ======================
%    DIAGRAMA DE BODE
% ======================
% Rango de frecuencias
f_min = 50;
f_max = 20000;

% Vector de evaluación
f = logspace(log10(f_min), log10(f_max), 2000);
w = 2 * pi * f;

% Respuesta en frecuencia
Hjw   = freqs(num, den, w);
magdB = 20*log10(abs(Hjw));
phase_deg = unwrap(angle(Hjw)) * 180/pi;
freqHz = f;  

figure('Units','normalized','Position',[0.12 0.12 0.65 0.6]);

% MAGNITUD
ax1 = gca;
set(ax1,'XScale','log'); hold(ax1,'on'); grid(ax1,'on');
semilogx(freqHz, magdB, 'b', 'LineWidth', 1.4);
xlabel('Frecuencia (Hz)');
ylabel('Magnitud (dB)');
title('Diagrama de Bode - Magnitud');

% Dibujar líneas horizontales de especificación
% Línea de -Ap (ondulación)
if exist('yline','file')==2
    yline(-Ap, '--g', sprintf('-A_p = %.2f dB', Ap), 'LabelHorizontalAlignment','left');
else
    semilogx([f_min, f_max], [-Ap, -Ap], '--g', 'LineWidth', 1);
end

if exist('yline','file')==2
    yline(-As, '--r', sprintf('-A_s = %d dB', As), 'LabelHorizontalAlignment','left');
else
    semilogx([f_min, f_max], [-As, -As], '--r', 'LineWidth', 1);
end

% Dibujar verticales en fs1, fp1, fp2, fs2
xline(fs1, ':k', sprintf('%d Hz', fs1)); 
xline(fp1, '--k', sprintf('%d Hz', fp1)); 
xline(fp2, '--k', sprintf('%d Hz', fp2)); 
xline(fs2, ':k', sprintf('%d Hz', fs2)); 

% Marcar y etiquetar puntos importantes
f_checks = [];
f_checks = [f_checks, fs1]; 
f_checks = [f_checks, fp1]; 
f_checks = [f_checks, fp2]; 
f_checks = [f_checks, fs2]; 

Hchk = freqs(num, den, 2*pi*f_checks);
magdB_chk = 20*log10(abs(Hchk));
plot(f_checks, magdB_chk, 'ko', 'MarkerFaceColor','k');
for k = 1:length(f_checks)
    text(f_checks(k)*1.05, magdB_chk(k), sprintf(' %g Hz: %.2f dB', f_checks(k), magdB_chk(k)), 'FontSize',8);
end

xlim([f_min f_max]);
ylim_auto = ylim;
ylim([max(-120, ylim_auto(1)), ylim_auto(2)]);

% Graficar cada sección SOS
nsecs = size(sos,1);
for sidx = 1:nsecs
    if sidx == 1
        bsec = 2*g*sos(sidx,1:3); asec = sos(sidx,4:6);
        Hsec = freqs(bsec, asec, w);
    else
        bsec = 1/2*sos(sidx,1:3); asec = sos(sidx,4:6);
        Hsec = freqs(bsec, asec, w);
    end
    semilogx(freqHz, 20*log10(abs(Hsec)), '--', 'LineWidth', 1);
end

legend('Filtro (magnitud)', 'Location','southwest');

%%
% =========================
%    PUNTOS IMPORTANTES
% =========================
f_check = [fs1 fp1 fp2 fs2];
w_check = 2 * pi * f_check;
Hchk = freqs(num, den, w_check);
magdB_chk = 20 * log10(abs(Hchk));
table(f_check', magdB_chk', 'VariableNames', {'f [Hz]','Mag [dB]'})

%%
% ========================================================
%    TOPOLOGÍA BICUADRÁTICA DE REALIMENTACIÓN POSITIVA
% ========================================================
% Definir Variables, siendo V1: Entrada Inversora de AO   (V-)
%                           V2: Tensión de Entrada        (Vin)
%                           V3: Tensión de Realimentación (Vout')
syms V1 V2 V3 Vx C1 C2 R1 R2 R3 k s

% Definir Sistema de Ecuaciones:
% --------------------------------------------------------------------------
% eq1 = (Vx - V2) / R1 + Vx * s * C1 + (Vx - V3) / R2 + (Vx - V1) * s * C2
% eq2 = V1 / R3 + (V1 - Vx) * s * C2
% --------------------------------------------------------------------------
% Hallar TFB (V2 = 0):
eq1 = Vx / R1 + Vx * s * C1 + (Vx - V3) / R2 + (Vx - V1) * s * C2;
eq2 = V1 / R3 + (V1 - Vx) * s * C2;

% Solución:
sol = solve(eq1, eq2, V1, Vx);
TFB = simplify(sol.V1 / V3);

% Hallar TFF (V3 = 0)
eq1 = (Vx - V2) / R1 + Vx * s * C1 + Vx / R2 + (Vx - V1) * s * C2;
eq2 = V1 / R3 + (V1 - Vx) * s * C2;

% Solución:
sol = solve(eq1, eq2, V1, Vx);
TFF = simplify(sol.V1 / V2);

% Hallar Tv (Función de Transferencia):
Tv = simplify(k * TFF / (1 - k * TFB));

%%
% ==========================
%    SINTESIS BLOQUE N°1
% ==========================

% ==============================
%    Función de Transferencia:
%            5252.26s
%   ----------------------------
%    s^2 + 3184.42s + 66320656
% ==============================

% De la FT se observa que:
wp1   = sqrt(66320656);
Qp1   = wp1 / 3184.42;

% Escalero de Impedancias:
C_b1     = 1;
R_b1     = sqrt(2) / wp1;
r2_r1_b1 = 3 - sqrt(2) / Qp1;

% Considerando 10^7:
C_b1 = C_b1 / 1E7;
R_b1 = R_b1 * 1E7;

% Por otro lado:
r1_b1 = 1E3;
r2_b1 = r2_r1_b1 * r1_b1;

% Resultados:
fprintf('            =================================\n');
fprintf('                 Componentes Bloque N°1      \n');
fprintf('            =================================\n');

fprintf('Capacitores:  C1 = C2 = %.0f [nF]\n', C_b1 * 1E9);
fprintf('Resistencias: R1 = R2 = R3 = %.2f [kΩ]\n', R_b1 / 1E3);
fprintf('Resistencias Realimentación: r1 = %.0f [kΩ]; r2 = %.3f [kΩ]\n', r1_b1 / 1E3, r2_b1 / 1E3);

%%
% ==========================
%    SINTESIS BLOQUE N°2
% ==========================

% ==============================
%    Función de Transferencia:
%            3126.49s
%   ----------------------------
%    s^2 + 1895.58s + 23500096
% ==============================

% De la FT se observa que:
wp2   = sqrt(23500096);
Qp2   = wp2 / 1895.58;

% Escalero de Impedancias:
C_b2     = 1;
R_b2     = sqrt(2) / wp2;
r2_r1_b2 = 3 - sqrt(2) / Qp2;

% Considerando 10^7:
C_b2 = C_b2 / 1E7;
R_b2 = R_b2 * 1E7;

% Por otro lado:
r1_b2 = 1E3;
r2_b2 = r2_r1_b2 * r1_b2;

% Resultados:
fprintf('            =================================\n');
fprintf('                 Componentes Bloque N°2      \n');
fprintf('            =================================\n');

fprintf('Capacitores:  C1 = C2 = %.0f [nF]\n', C_b2 * 1E9);
fprintf('Resistencias: R1 = R2 = R3 = %.2f [kΩ]\n', R_b2 / 1E3);
fprintf('Resistencias Realimentación: r1 = %.0f [kΩ]; r2 = %.3f [kΩ]\n', r1_b2 / 1E3, r2_b2 / 1E3);

%%
% ===================
%    SENSIBILIDAD 
% ===================
% De la FT se tiene que:
wp  = simplify(sqrt((R1 + R2) / (C1 * C2 * R1 * R2 * R3)));
bwp = simplify((C1 * R1 * R2 + C2 * R1 * R2 + C2 * R1 * R3 + C2 * R2 * R3 - C2 * R1 * R3 * k) / (C1 * C2 * R1 * R2 * R3));

% Sensibilidad de cada Componente:
S_wp_R1 = simplify(diff(wp, R1) * R1 / wp);
S_wp_R2 = simplify(diff(wp, R2) * R2 / wp);
S_wp_R3 = simplify(diff(wp, R3) * R3 / wp);
S_wp_C1 = simplify(diff(wp, C1) * C1 / wp);
S_wp_C2 = simplify(diff(wp, C1) * C1 / wp);

S_bwp_R1 = simplify(diff(bwp, R1) * R1 / bwp);
S_bwp_R2 = simplify(diff(bwp, R2) * R2 / bwp);
S_bwp_R3 = simplify(diff(bwp, R3) * R3 / bwp);
S_bwp_C1 = simplify(diff(bwp, C1) * C1 / bwp);
S_bwp_C2 = simplify(diff(bwp, C2) * C2 / bwp);

% Recordar que: R1 = R2 = R3 y C1 = C2:
syms R C

% Resultados:
fprintf(' ==============================\n');
fprintf('         Sensibilidad        \n');
fprintf(' ==============================\n');

fprintf('wp:\n');
fprintf('       S_wp_R1 = %.2f\n', double(subs(S_wp_R1, [R1, R2, R3, C1, C2], [R, R, R, C, C])));
fprintf('       S_wp_R2 = %.2f\n', double(subs(S_wp_R2, [R1, R2, R3, C1, C2], [R, R, R, C, C])));
fprintf('       S_wp_R3 = %.1f\n', double(subs(S_wp_R3, [R1, R2, R3, C1, C2], [R, R, R, C, C])));
fprintf('       S_wp_C1 = %.1f\n', double(subs(S_wp_C1, [R1, R2, R3, C1, C2], [R, R, R, C, C])));
fprintf('       S_wp_C2 = %.1f\n', double(subs(S_wp_C2, [R1, R2, R3, C1, C2], [R, R, R, C, C])));
fprintf(' ------------------------------\n');
fprintf('bwp:\n');
fprintf('      S_bwp_R1 = %.2f\n', double(subs(S_bwp_R1, [R1, R2, R3, C1, C2, k], [R, R, R, C, C, r2_r1_b1 + 1])));
fprintf('      S_bwp_R2 = %.2f\n', double(subs(S_bwp_R2, [R1, R2, R3, C1, C2, k], [R, R, R, C, C, r2_r1_b1 + 1])));
fprintf('      S_bwp_R3 = %.2f\n', double(subs(S_bwp_R3, [R1, R2, R3, C1, C2, k], [R, R, R, C, C, r2_r1_b1 + 1])));
fprintf('      S_bwp_C1 = %.2f\n', double(subs(S_bwp_C1, [R1, R2, R3, C1, C2, k], [R, R, R, C, C, r2_r1_b1 + 1])));
fprintf('      S_bwp_C2 = %.2f\n', double(subs(S_bwp_C2, [R1, R2, R3, C1, C2, k], [R, R, R, C, C, r2_r1_b1 + 1])));
