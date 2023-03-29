close all;

d_s = 1.5e-3;
D_s = d_s + 4.75e-3;
N   = 25;
R   = D_s / d_s;
W_c = (4 * R - 1) / (4 * R - 4) + 0.615 / R;

tau_scr = 60e6;
tau_fcr = 320e6;
G_M     = 10.39e9;
% G_M     = 14.39e9;
G_A     = 16.05e9;
% G_A     = 18.75e9;
gamma_L = 0.02;
delta_L = gamma_L * pi * N * D_s^2 / d_s;
C1 = 8 * W_c * D_s / (pi * d_s^3);
C2 = d_s / (pi * N * D_s^2);
K = C2 / C1;%d_s^4 / (8 * N * W_c * D_s^3);
Theta = 1.504e6;

A_s = 29.23;
% A_s = 40.83;
A_f = 51.55;
% A_f = 44.55;
C_M = 11.51e6;
C_A = 36.05e6;
% C_A = 39.55e6;
a_A = pi / (A_f - A_s);

% calculate force and xi_s during pre-stretch
delta = (0 : 5 : 55) * 1e-3;
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
T = 22 : 0.2 : 60;
xi        = zeros(length(delta), length(T));
xi_S_heat = zeros(length(delta), length(T));
x_heat    = zeros(length(delta), length(T));
% for i = 1 : length(delta)
%     x_heat(i, 1) = x_pre(i);
%     for j = 2 : length(T)
%         if T(j) < A_s + x_heat(i, j - 1) * C1 / C_A
%             x_heat(i, j) = fsolve( ...
%             @(F) C1 * (F - x_pre(i)), x_heat(i, j-1));
%         elseif T(j) < A_f + x_heat(i, j - 1) * C1 / C_A
%             x_heat(i, j) = fsolve( ...
%                 @(F) C1 * (F - x_pre(i)) - C2 * ( ...
%                     (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * delta(i) - G_M * delta(i) ...
%                 ) + C2 * delta_L * ( ...
%                     (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * ...
%                         (xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1))) - ...
%                     G_M * xi_S_pre(i) ...
%                 ), ...
%                 x_heat(i, j-1));
%         else
%             x_heat(i, j) = fsolve( ...
%             @(F) C1 * (F - x_pre(i)) - C2 * ( ...
%                     G_A * delta(i) - G_M * delta(i) ...
%                 ) + C2 * delta_L * ( ...
%                     G_A * 0 - G_M * xi_S_pre(i) ...
%                 ), ...
%                 x_heat(i, j-1));
%         end
%     end
% end

for i = 1 : length(delta)
    x_heat(i, 1) = x_pre(i);
    for j = 2 : length(T)
        if T(j) < A_s + x_heat(i, j - 1) * C1 / C_A
            x_heat(i, j) = fsolve( ...
            @(F) C1 * (F - x_pre(i)) - Theta * (T(j) - T(1)), x_heat(i, j-1));
        elseif T(j) < A_f + x_heat(i, j - 1) * C1 / C_A
            x_heat(i, j) = fsolve( ...
                @(F) C1 * (F - x_pre(i)) - C2 * ( ...
                    (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * delta(i) - G_M * delta(i) ...
                ) + C2 * delta_L * ( ...
                    (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1)) * ...
                        (xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - F*C1/C_A)) + 1))) - ...
                    G_M * xi_S_pre(i) ...
                ) - Theta * (T(j) - T(1)), ...
                x_heat(i, j-1));
        else
            x_heat(i, j) = fsolve( ...
            @(F) C1 * (F - x_pre(i)) - C2 * ( ...
                    G_A * delta(i) - G_M * delta(i) ...
                ) + C2 * delta_L * ( ...
                    G_A * 0 - G_M * xi_S_pre(i) ...
                ) - Theta * (T(j) - T(1)), ...
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
% figure('Name', 'block force');
% hold on;
% grid on;
% box on;
% for i = length(delta) : -1 : 6
%     plot(T, (x_heat(i, :) - x_pre(i)) * 5.2, '-', 'LineWidth', 1);%, 'Color', color_map(i - 3, :) / 255
% end
% legend('50mm', '45mm', '40mm', '35mm', '30mm', 'Location', 'northwest');
% hold off;
% xlabel('Temperature ({}^{\circ}C)', 'Fontsize', 8);
% ylabel('Torque (Nmm)', 'Fontsize', 8);
% xlim([20 65]);
% ylim([0 200]);
% % xticks(0.1 : 0.1 : 0.5);
% % yticks([-90 -60 -30 0]);
% set(gca, 'FontName', 'Times New Roman','Fontsize', 8);

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 1.5], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 1.5]);
% saveas(gcf, '..\..\paper\v1.0.0\figures\blocking_force_vs_temperature.pdf');
% filename = '..\..\paper\v1.0.0\figures\blocking_force_vs_temperature.tiff';
% % print(gcf, filename, '-dtiffn', '-r600');

%------------------------------------------------------------------------------------------%
% plot stall torque vs temperature at different predislacement
close all;
% clear;

% directory check
% data_file_path = checkDirectory();

datafilelist = { ...
'2021-01-08 11-02-08_smpl. freq=500Hz,stall torque,L0=20mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 11-44-39_smpl. freq=500Hz,stall torque,L0=25mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 12-17-00_smpl. freq=500Hz,stall torque,L0=30mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 13-39-16_smpl. freq=500Hz,stall torque,L0=35mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 14-17-16_smpl. freq=500Hz,stall torque,L0=40mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 14-53-32_smpl. freq=500Hz,stall torque,L0=45mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 15-29-32_smpl. freq=500Hz,stall torque,L0=50mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-08 16-01-05_smpl. freq=500Hz,stall torque,L0=55mm,passive cooling,open loop,0.25A__data.mat' ...
};
datafilelist_F = { ...
'2021-01-08 11-02-08_smpl. freq=500Hz,stall torque,L0=20mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 11-44-39_smpl. freq=500Hz,stall torque,L0=25mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 12-17-00_smpl. freq=500Hz,stall torque,L0=30mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 13-39-16_smpl. freq=500Hz,stall torque,L0=35mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 14-17-16_smpl. freq=500Hz,stall torque,L0=40mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 14-53-32_smpl. freq=500Hz,stall torque,L0=45mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 15-29-32_smpl. freq=500Hz,stall torque,L0=50mm,passive cooling,open loop,0.25A__F_and_T.mat', ...
'2021-01-08 16-01-05_smpl. freq=500Hz,stall torque,L0=55mm,passive cooling,open loop,0.25A__F_and_T.mat' ...
};
% datafilelist   = strcat(data_file_path, datafilelist);
% datafilelist_F = strcat(data_file_path, datafilelist_F);
sampling_frequency = 500;
radius_transducer  = 40;
predislacement     = (20 : 5 : 55)';
torque_cell        = cell(1, length(datafilelist_F));
T_cell             = cell(1, length(datafilelist_F));
end_point          = [4086, 3527, 4490, 5827, 5743, 6223, 5822, 4470];

% extract force and temperature (sampling frequency is 1/10 of temeprature data)
% for i = 1 : length(datafilelist_F)
%     load(datafilelist_F{i});
%     % torque_cell{i} = ...
%     %     sqrt( ...
%     %         F_and_T(2, 101 : end_point(i)).^2 + ...
%     %         F_and_T(3, 101 : end_point(i)).^2 + ...
%     %         F_and_T(4, 101 : end_point(i)).^2 ...
%     %     );
%     torque_cell{i} = F_and_T(4, 101 : end_point(i));
%     torque_cell{i} = torque_cell{i} - torque_cell{i}(1);
%     torque_cell{i} = torque_cell{i} * radius_transducer;

%     load(datafilelist{i});
%     T_cell{i} = data(2, 1 : 100 : (length(torque_cell{i})) * 100);
% end


% plot torque vs temperature
figure('Name', 'stall torque');
hold on;
% plot(T_cell{5}, abs(torque_cell{5}), 'LineWidth', 1, 'Color', [  0  114  189] / 255);
plot(T, (x_heat(9, :) - x_pre(9)) * 5.2, '--', 'LineWidth', 1, 'Color', [  0  114  189] / 255);

% plot(T_cell{4}, abs(torque_cell{4}), 'LineWidth', 1, 'Color', [237  177   32] / 255);
plot(T, (x_heat(8, :) - x_pre(8)) * 5.2, '--', 'LineWidth', 1, 'Color', [237  177   32] / 255);

% plot(T_cell{2}, abs(torque_cell{2}), 'LineWidth', 1, 'Color', [126   47  142] / 255);
plot(T, (x_heat(7, :) - x_pre(7)) * 5.2, '--', 'LineWidth', 1, 'Color', [126   47  142] / 255);
hold off;
% leg = legend('40 mm (Expt.)', '40 mm (Model)', '35 mm (Expt.)', '35 mm (Model)', '30 mm (Expt.)', '30 mm (Model)', 'NumColumns', 3, 'Location', 'northout');
% leg.ItemTokenSize = [20, 20];
xlabel('$T_{s}$\,\,($^{\circ}$C)', 'interpreter', 'latex', 'Fontsize', 8);
ylabel('$T_{st}$\,\,(Nmm)', 'interpreter', 'latex', 'Fontsize', 8);
xlim([20 60]);
ylim([0 200]);
xticks(20 : 10 : 60);
yticks(0 : 50 : 200);
grid on;
box on;
set(gca, 'FontName', 'Times New Roman','Fontsize', 8);

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 1.75, 1.3], 'PaperUnits', 'Inches', 'PaperSize', [1.75, 1.3]);
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 1.3], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 1.3]);
% saveas(gcf, '..\..\paper\v1.0.0\figures\stall_torque_vs_temperature.pdf');
filename = '..\..\paper\v1.0.0\figures\stall_torque_vs_temperature.tiff';
% print(gcf, filename, '-dtiffn', '-r600');
