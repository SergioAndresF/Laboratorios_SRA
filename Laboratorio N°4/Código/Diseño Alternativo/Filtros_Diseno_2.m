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
syms V1 V2 V3 Vx C1 C2 R1 R2 k s

% Definir Sistema de Ecuaciones:
% --------------------------------------------------------------------------
% eq1 = (Vx - V2) / R1 + (Vx - V3) * s * C1 + (Vx - V1) / R2
% eq2 = V1 * s * C2 + (V1 - Vx) / R2
% --------------------------------------------------------------------------
% Hallar TFB (V2 = 0):
eq1 = Vx / R1 + (Vx - V3) * s * C1 + (Vx - V1) / R2;
eq2 = V1 * s * C2 + (V1 - Vx) / R2;

% Solución:
sol    = solve(eq1, eq2, V1, Vx);
TFB_b1 = simplify(sol.V1 / V3);

% Hallar TFF (V3 = 0)
eq1 = (Vx - V2) / R1 + Vx * s * C1 + (Vx - V1) / R2;
eq2 = V1 * s * C2 + (V1 - Vx) / R2;

% Solución:
sol    = solve(eq1, eq2, V1, Vx);
TFF_b1 = simplify(sol.V1 / V2);

% Hallar Tv (Función de Transferencia):
Tv_b1 = simplify(k * TFF_b1 / (1 - k * TFB_b1));

%%
% ==========================
%    SINTESIS BLOQUE N°1
% ==========================

% ==============================
%    Función de Transferencia:
%            3.284e07
%   ---------------------------- 
%     s^2 + 3184 s + 6.632e07
% ==============================

% De la FT se observa que:
wp1   = sqrt(6.632E07);
Qp1   = wp1 / 3184;

% Escalero de Impedancias:
C_b1     = 1;
R1_b1    = 1 / 3184;
R2_b1    = 3184 / 6.632e07;

% Considerando 10^7:
C_b1  = C_b1 / 1E7;
R1_b1 = R1_b1 * 1E7;
R2_b1 = R2_b1 * 1E7;
r_b1  = 1E3;

% Resultados:
fprintf('            =================================\n');
fprintf('                 Componentes Bloque N°1      \n');
fprintf('            =================================\n');

fprintf('Capacitores:  C1 = C2 = %.0f [nF]\n', C_b1 * 1E9);
fprintf('Resistencias: R1 = %.2f [kΩ]; R2 = %.0f [Ω]\n', R1_b1 / 1E3, R2_b1);
fprintf('Resistencias Realimentación: r1 = r2 = %.f [kΩ]\n', r_b1 / 1E3);

%%
% ===================
%    SENSIBILIDAD 
% ===================
% De la FT se tiene que:
wp_b1  = simplify(sqrt(1 / (C1 * C2 * R1 * R2)));
bwp_b1 = simplify((C1 * R1 + C2 * R1 + C2 * R2 - C1 * R1 * k) / (C1 * C2 * R1 * R2));

% Sensibilidad de cada Componente:
S_wp_R1 = simplify(diff(wp_b1, R1) * R1 / wp_b1);
S_wp_R2 = simplify(diff(wp_b1, R2) * R2 / wp_b1);
S_wp_C1 = simplify(diff(wp_b1, C1) * C1 / wp_b1);
S_wp_C2 = simplify(diff(wp_b1, C1) * C1 / wp_b1);

S_bwp_R1 = simplify(diff(bwp_b1, R1) * R1 / bwp_b1);
S_bwp_R2 = simplify(diff(bwp_b1, R2) * R2 / bwp_b1);
S_bwp_C1 = simplify(diff(bwp_b1, C1) * C1 / bwp_b1);
S_bwp_C2 = simplify(diff(bwp_b1, C2) * C2 / bwp_b1);

% Recordar que: R1 = R2 = R3 y C1 = C2:
syms R C

% Resultados:
fprintf(' ==============================\n');
fprintf('         Sensibilidad        \n');
fprintf(' ==============================\n');

fprintf('wp:\n');
fprintf('       S_wp_R1 = %.1f\n', double(subs(S_wp_R1, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_R2 = %.1f\n', double(subs(S_wp_R2, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_C1 = %.1f\n', double(subs(S_wp_C1, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_C2 = %.1f\n', double(subs(S_wp_C2, [R1, R2, C1, C2], [R, R, C, C])));
fprintf(' ------------------------------\n');
fprintf('bwp:\n');
fprintf('      S_bwp_R1 = %.1f\n', double(subs(S_bwp_R1, [R1, R2, C1, C2, k], [R, R, C, C, 2])));
fprintf('      S_bwp_R2 = %.1f\n', double(subs(S_bwp_R2, [R1, R2, C1, C2, k], [R, R, C, C, 2])));
fprintf('      S_bwp_C1 = %.1f\n', double(subs(S_bwp_C1, [R1, R2, C1, C2, k], [R, R, C, C, 2])));
fprintf('      S_bwp_C2 =  %.1f\n', double(subs(S_bwp_C2, [R1, R2, C1, C2, k], [R, R, C, C, 2])));

%%
% ========================================================
%    TOPOLOGÍA BICUADRÁTICA DE REALIMENTACIÓN POSITIVA
% ========================================================
% Definir Variables, siendo V1: Entrada Inversora de AO   (V-)
%                           V2: Tensión de Entrada        (Vin)
%                           V3: Tensión de Realimentación (Vout')
syms V1 V2 V3 Vx C1 C2 R1 R2 k s

% Definir Sistema de Ecuaciones:
% --------------------------------------------------------------------------
% eq1 = (Vx - V2) * s * C1 + (Vx - V3) / R1 + (Vx - V1) * s * C2
% eq2 = V1 / R2 + (V1 - Vx) * s * C2
% --------------------------------------------------------------------------
% Hallar TFB (V2 = 0):
eq1 = Vx * s * C1 + (Vx - V3) / R1 + (Vx - V1) * s * C2;
eq2 = V1 / R2 + (V1 - Vx) * s * C2;

% Solución:
sol    = solve(eq1, eq2, V1, Vx);
TFB_b2 = simplify(sol.V1 / V3);

% Hallar TFF (V3 = 0)
eq1 = (Vx - V2) * s * C1 + Vx / R1 + (Vx - V1) * s * C2;
eq2 = V1 / R2 + (V1 - Vx) * s * C2;

% Solución:
sol    = solve(eq1, eq2, V1, Vx);
TFF_b2 = simplify(sol.V1 / V2);

% Hallar Tv (Función de Transferencia):
Tv_b2 = simplify(k * TFF_b2 / (1 - k * TFB_b2));

%%
% ==========================
%    SINTESIS BLOQUE N°2
% ==========================

% ==============================
%    Función de Transferencia:
%            0.5 s^2
%   ----------------------------
%    s^2 + 1896 s + 2.35e07
% ==============================

% De la FT se observa que:
wp2   = sqrt(2.35E7);
Qp2   = wp2 / 1896;

% Escalero de Impedancias:
C_b2     = 1;
R_b2     = 1 / wp2;
r2_r1_b2 = 2 - wp2 / Qp2 * R_b2;

% Considerando 10^7:
C_b2  = C_b2 / 1E7;
R_b2  = R_b2 * 1E7;
r1_b2 = 1E3;
r2_b2 = r2_r1_b2 * r1_b2;

% Resultados:
fprintf('            =================================\n');
fprintf('                 Componentes Bloque N°2      \n');
fprintf('            =================================\n');

fprintf('Capacitores:  C1 = C2 = %.0f [nF]\n', C_b2 * 1E9);
fprintf('Resistencias: R1 = R2 = %.2f [kΩ]\n', R_b2 / 1E3);
fprintf('Resistencias Realimentación: r1 = %.0f [kΩ]; r2 = %.3f [kΩ]\n',r1_b2 / 1E3, r2_b2 / 1E3);

%%
% ===================
%    SENSIBILIDAD 
% ===================
% De la FT se tiene que:
wp_b2  = simplify(sqrt(1 / (C1 * C2 * R1 * R2)));
bwp_b2 = simplify((C1 * R1 + C2 * R1 + C2 * R2 - C2 * R2 * k) / (C1 * C2 * R1 * R2));

% Sensibilidad de cada Componente:
S_wp_R1 = simplify(diff(wp_b2, R1) * R1 / wp_b2);
S_wp_R2 = simplify(diff(wp_b2, R2) * R2 / wp_b2);
S_wp_C1 = simplify(diff(wp_b2, C1) * C1 / wp_b2);
S_wp_C2 = simplify(diff(wp_b2, C2) * C2 / wp_b2);

S_bwp_R1 = simplify(diff(bwp_b2, R1) * R1 / bwp_b2);
S_bwp_R2 = simplify(diff(bwp_b2, R2) * R2 / bwp_b2);
S_bwp_C1 = simplify(diff(bwp_b2, C1) * C1 / bwp_b2);
S_bwp_C2 = simplify(diff(bwp_b2, C2) * C2 / bwp_b2);

% Recordar que: R1 = R2 = R3 y C1 = C2:
syms R C

% Resultados:
fprintf(' ==============================\n');
fprintf('         Sensibilidad        \n');
fprintf(' ==============================\n');

fprintf('wp:\n');
fprintf('       S_wp_R1 = %.1f\n', double(subs(S_wp_R1, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_R2 = %.1f\n', double(subs(S_wp_R2, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_C1 = %.1f\n', double(subs(S_wp_C1, [R1, R2, C1, C2], [R, R, C, C])));
fprintf('       S_wp_C2 = %.1f\n', double(subs(S_wp_C2, [R1, R2, C1, C2], [R, R, C, C])));
fprintf(' ------------------------------\n');
fprintf('bwp:\n');
fprintf('      S_bwp_R1 =  %.2f\n', double(subs(S_bwp_R1, [R1, R2, C1, C2, k], [R, R, C, C, r2_r1_b2 + 1])));
fprintf('      S_bwp_R2 = %.2f\n', double(subs(S_bwp_R2, [R1, R2, C1, C2, k], [R, R, C, C, r2_r1_b2 + 1])));
fprintf('      S_bwp_C1 =  %.2f\n', double(subs(S_bwp_C1, [R1, R2, C1, C2, k], [R, R, C, C, r2_r1_b2 + 1])));
fprintf('      S_bwp_C2 = %.2f\n', double(subs(S_bwp_C2, [R1, R2, C1, C2, k], [R, R, C, C, r2_r1_b2 + 1])));
