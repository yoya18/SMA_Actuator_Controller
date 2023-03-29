close all;

d_s = 1.5e-3;
D_s = d_s + 4.75e-3;
N   = 25;
R   = D_s / d_s;
W_c = (4 * R - 1) / (4 * R - 4) + 0.615 / R;

tau_scr = 60e6;
tau_fcr = 320e6;
G_M     = 14.39e9;
G_A     = 18.75e9;
gamma_L = 0.02;
delta_L = gamma_L * pi * N * D_s^2 / d_s;
C1 = 8 * W_c * D_s / (pi * d_s^3);
C2 = d_s / (pi * N * D_s^2);
K = C2 / C1;%d_s^4 / (8 * N * W_c * D_s^3);

A_s = 40.83;
A_f = 44.55;
C_M = 11.51e6;
C_A = 39.55e6;
a_A = pi / (A_f - A_s);

% calculate force and xi_s during pre-stretch
delta = (0 : 5 : 50) * 1e-3;
xi_S_pre  = zeros(1, length(delta));
x_pre = zeros(1, length(delta));
for i = 1 : length(delta)
    if i > 1
        if x_pre(i - 1) < tau_scr / C1
            xi_S_pre(i) = 0;
            x_pre(i) = fsolve( ...
            @(F) C1 * F - C2 * G_M * delta(i) + ...
                C2 * delta_L * G_M * (xi_S_pre(i)), ...
            40);
        elseif x_pre(i - 1) < tau_fcr / C1
            x_pre(i) = fsolve( ...
                @(F) C1 * F - C2 * G_M * delta(i) + ...
                    C2 * delta_L * G_M * (1/2 * cos(pi / (tau_scr - tau_fcr) * (F*C1 - tau_fcr)) + 1/2), ...
                40);
            xi_S_pre(i) = 1/2 * cos(pi / (tau_scr - tau_fcr) * (x_pre(i)*C1 - tau_fcr)) + 1/2;
        else
            x_pre(i) = fsolve( ...
            @(F) C1 * F - C2 * G_M * delta(i) + ...
                C2 * delta_L * G_M * (1), ...
            40);
            xi_S_pre(i) = 1;
        end
    end
end

% calculate force of active spring
T = 20 : 0.2 : 60;
xi        = zeros(length(delta), length(T));
xi_S_heat = zeros(length(delta), length(T));
x_heat    = zeros(length(delta), length(T));
temp = zeros(length(delta), length(T));
for i = 1 : length(delta)
    x_heat(i, 1) = x_pre(i);
    for j = 2 : length(T)
        if T(j) < A_s + x_heat(i, j - 1) * C1 / C_A
            x_heat(i, j) = fsolve( ...
            @(F) C1 * (F - x_pre(i)), x_heat(i, j-1));
        elseif T(j) < A_f + x_heat(i, j - 1) * C1 / C_A
            x_heat(i, j) = fsolve( ...
                @(F) C1 * (F - x_pre(i)) - C2 * ( ...
                    (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * delta(i) - G_M * delta(i) ...
                ) + C2 * delta_L * ( ...
                    (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * ...
                        (xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1))) - ...
                    G_M * xi_S_pre(i) ...
                ), ...
                x_heat(i, j-1));
        else
            x_heat(i, j) = fsolve( ...
            @(F) C1 * (F - x_pre(i)) - C2 * ( ...
                    G_A * delta(i) - G_M * delta(i) ...
                ) + C2 * delta_L * ( ...
                    G_A * 0 - G_M * xi_S_pre(i) ...
                ), ...
                x_heat(i, j-1));
        end
    end
end

% figure('Name', 'force vs. dis');
% hold on;
% grid on;
% plot(delta, x_pre);
% plot([0, delta(end)], [tau_scr / C1, tau_scr / C1]);
% plot([0, delta(end)], [tau_fcr / C1, tau_fcr / C1]);
% % plot(delta * C2, x_pre * C1);
% % plot([0, delta(end) * C2], [tau_scr, tau_scr]);
% % plot([0, delta(end) * C2], [tau_fcr, tau_fcr]);
% hold off;

% figure('Name', 'stress vs. strain');
% hold on;
% grid on;
% plot(delta * C2, x_pre * C1 * 1e-6);
% plot([0, delta(end) * C2], [tau_scr * 1e-6, tau_scr * 1e-6]);
% plot([0, delta(end) * C2], [tau_fcr * 1e-6, tau_fcr * 1e-6]);
% hold off;

color_map = [ ...
    217   83   25; ...
      0  114  189; ...
    237  177   32; ...
    126   47  142; ...
    119  177   48
];
figure('Name', 'block force');
hold on;
grid on;
box on;
for i = length(delta) : -1 : 7
    plot(T, (x_heat(i, :) - x_pre(i)) * 5.2, '-', 'LineWidth', 1, 'Color', color_map(i - 6, :) / 255);
end
legend('50mm', '45mm', '40mm', '35mm', '30mm', 'Location', 'northwest');
hold off;
xlabel('Temperature ({}^{\circ}C)', 'Fontsize', 8);
ylabel('Torque (Nmm)', 'Fontsize', 8);
xlim([40 60]);
ylim([0 250]);
% xticks(0.1 : 0.1 : 0.5);
% yticks([-90 -60 -30 0]);
set(gca, 'FontName', 'Times New Roman','Fontsize', 8);

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 1.5], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 1.5]);
saveas(gcf, '..\..\paper\v1.0.0\figures\blocking_force_vs_temperature.pdf');
filename = '..\..\paper\v1.0.0\figures\blocking_force_vs_temperature.tiff';
% print(gcf, filename, '-dtiffn', '-r600');