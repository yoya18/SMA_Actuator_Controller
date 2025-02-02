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
G_A     = 11.25e9;
% G_A     = 18.75e9;
gamma_L = 0.02;
delta_L = gamma_L * pi * N * D_s^2 / d_s;
C1 = 8 * W_c * D_s / (pi * d_s^3);
C2 = d_s / (pi * N * D_s^2);
K = C2 / C1;%d_s^4 / (8 * N * W_c * D_s^3);
Theta = 0.704e6;

A_s = 34.23;
% A_s = 40.83;
A_f = 55.55;
% A_f = 44.55;
C_M = 11.51e6;
C_A = 30.05e6;
% C_A = 39.55e6;
a_A = pi / (A_f - A_s);

% calculate force and xi_s during pre-stretch
delta = (0 : 5 : 55) * 1e-3;
xi_S_pre = zeros(1, length(delta));
x_pre = zeros(1, length(delta));
for i = 1 : length(delta)
    % if delta(i) == 30
    %     G_A = 11.25e9;
    % elseif delta(i) == 35
    %     G_A = 14.25e9;
    % elseif delta(i) == 40
    %     G_A = 18.25e9;
    % end
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
T = 30 : 0.2 : 60;
x_heat = zeros(length(delta), length(T));
l      = zeros(length(delta), length(T));
xi     = zeros(length(delta), length(T));
for i = 1 : length(delta)
    if delta(i) == 30e-3
        G_A = 11.25e9;
    elseif delta(i) == 35e-3
        G_A = 12.25e9;
    elseif delta(i) == 40e-3
        G_A = 12.45e9;
    end
    x_heat(i, 1) = x_pre(i);
    l(i, 1) = delta(i);
    for j = 2 : length(T)
        if T(j) < A_s + x_heat(i, j - 1) * C1 / C_A
            % x_heat(i, j) = x_pre(i) + Theta * (T(j) - T(1));
            % l(i, j)      = delta(i);
            if x_heat(i, j - 1) * C1 < tau_scr
                fucn = @(Fl) [ ...
                    Fl(1) - K * G_M * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            xi_S_pre(i) ...
                        ) ...
                    )
                ];
            elseif x_heat(i, j - 1) * C1 < tau_fcr
                fucn = @(Fl) [ ...
                    Fl(1) - K * G_M * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            1/2 * cos(pi / (tau_scr - tau_fcr) * (Fl(1)*C1 - tau_fcr)) + 1/2 ...
                        ) ...
                    )
                ];
            else
                fucn = @(Fl) [ ...
                    Fl(1) - K * G_M * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            1 ...
                        ) ...
                    )
                ];
            end
            S = fsolve(fucn, [x_heat(i, j-1); l(i, j-1)]);
            x_heat(i, j) = S(1);
            l(i, j)      = S(2);
        elseif T(j) < A_f + x_heat(i, j - 1) * C1 / C_A
            % if x_heat(i, j - 1) * C1 < tau_scr
            %     fucn = @(Fl) [ ...
            %         Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
            %             Fl(2) - delta_L * ( ...
            %                 xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
            %             ) ...
            %         );
            %         Fl(1) - K * G_M * ( ...
            %             2 * delta(i) - Fl(2) - delta_L * ( ...
            %                 xi_S_pre(i) ...
            %             ) ...
            %         )
            %     ];
            % elseif x_heat(i, j - 1) * C1 < tau_fcr
            %     fucn = @(Fl) [ ...
            %         Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
            %             Fl(2) - delta_L * ( ...
            %                 xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
            %             ) ...
            %         );
            %         Fl(1) - K * G_M * ( ...
            %             2 * delta(i) - Fl(2) - delta_L * ( ...
            %                 (1 - xi_S_pre(i))/2 * cos(pi / (tau_scr - tau_fcr) * (Fl(1)*C1 - tau_fcr)) + (1 + xi_S_pre(i))/2 ...
            %             ) ...
            %         )
            %     ];
            % else
            %     fucn = @(Fl) [ ...
            %         Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
            %             Fl(2) - delta_L * ( ...
            %                 xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
            %             ) ...
            %         );
            %         Fl(1) - K * G_M * ( ...
            %             2 * delta(i) - Fl(2) - delta_L * ( ...
            %                 1 ...
            %             ) ...
            %         )
            %     ];
            % end

            if x_heat(i, j - 1) * C1 < tau_scr
                fucn = @(Fl) [ ...
                    Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            xi_S_pre(i) ...
                        ) ...
                    )
                ];
            elseif x_heat(i, j - 1) * C1 < tau_fcr
                fucn = @(Fl) [ ...
                    Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            1/2 * cos(pi / (tau_scr - tau_fcr) * (Fl(1)*C1 - tau_fcr)) + 1/2 ...
                        ) ...
                    )
                ];
            else
                fucn = @(Fl) [ ...
                    Fl(1) - K * (G_A + (G_M - G_A) * 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) * ( ...
                        Fl(2) - delta_L * ( ...
                            xi_S_pre(i) - xi_S_pre(i)/1 * (1 - 1/2 * (cos(a_A * (T(j) - A_s - Fl(1)*C1/C_A)) + 1)) ...
                        ) ...
                    ) - Theta / C1 * (T(j) - T(1));
                    Fl(1) - K * G_M * ( ...
                        2 * delta(i) - Fl(2) - delta_L * ( ...
                            1 ...
                        ) ...
                    )
                ];
            end

            S = fsolve(fucn, [x_heat(i, j-1); l(i, j-1)]);
            x_heat(i, j) = S(1);
            l(i, j)      = S(2);
            % xi(i, j) = 1/2 * (cos(a_A * (T(j) - A_s - x_heat(i, j)*C1/C_A)) + 1);
        else
            x_heat(i, j) = x_heat(i, j-1);
            l(i, j)      = l(i, j-1);
        end
    end
end

% color_map = [ ...
%     217   83   25; ...
%       0  114  189; ...
%     237  177   32; ...
%     126   47  142; ...
%     119  177   48
% ];
% figure('Name', 'stroke');
% hold on;
% grid on;
% box on;
% for i = length(delta) : -1 : 7
%     plot(T, (delta(i) - l(i, :)) * 1e3 / 5.2 / pi * 180, '-', 'LineWidth', 1, 'Color', color_map(i - 6, :) / 255);
% end
% legend('50mm', '45mm', '40mm', '35mm', '30mm', 'Location', 'northwest');
% hold off;
% xlabel('Temperature ({}^{\circ}C)', 'Fontsize', 8);
% ylabel('Stroke (degree)', 'Fontsize', 8);
% xlim([40 60]);
% ylim([0 260]);
% % xticks(0.1 : 0.1 : 0.5);
% yticks(0 : 50 : 250);
% set(gca, 'FontName', 'Times New Roman','Fontsize', 8);

% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3.5, 1.5], 'PaperUnits', 'Inches', 'PaperSize', [3.5, 1.5]);
% % saveas(gcf, '..\..\paper\v1.0.0\figures\stroke_vs_temperature.pdf');
% filename = '..\..\paper\v1.0.0\figures\stroke_vs_temperature.tiff';
% % print(gcf, filename, '-dtiffn', '-r600');


%------------------------------------------------------------------------------------------%
% plot stall torque vs temperature at different predislacement
close all;
% clear;

% directory check
data_file_path = checkDirectory();

datafilelist = { ...
'2021-01-10 15-06-59_smpl. freq=500Hz,stroke,L0=20mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 15-38-26_smpl. freq=500Hz,stroke,L0=25mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 16-08-31_smpl. freq=500Hz,stroke,L0=30mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 16-36-11_smpl. freq=500Hz,stroke,L0=35mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 17-06-05_smpl. freq=500Hz,stroke,L0=40mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 17-32-00_smpl. freq=500Hz,stroke,L0=45mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 17-55-09_smpl. freq=500Hz,stroke,L0=50mm,passive cooling,open loop,0.25A__data.mat', ...
'2021-01-10 18-17-51_smpl. freq=500Hz,stroke,L0=55mm,passive cooling,open loop,0.25A__data.mat' ...
};
datafilelist   = strcat(data_file_path, datafilelist);
sampling_frequency = 500;
radius_transducer  = 40;
predislacement     = (20 : 5 : 55)';
angle_cell         = cell(1, length(datafilelist));
T_cell             = cell(1, length(datafilelist));
end_point          = [338200, 358300, 373200, 424000, 425400, 402900, 288100, 339900];

% extract force and temperature (sampling frequency is 1/10 of temeprature data)
for i = 1 : length(datafilelist)
    load(datafilelist{i});
    angle_cell{i} = data(10, 200 : 100 : end_point(i));
    T_cell{i} = data(2, 200 : 100 : end_point(i));
end

% [  0  114  189] blue
% [217   83   25] red
% colorMAP_R = linspace(  0, 217, length(current)).';
% colorMAP_G = linspace(114,  83, length(current)).';
% colorMAP_B = linspace(189,  25, length(current)).';
% colorMAP   = [colorMAP_R, colorMAP_G, colorMAP_B];

% plot torque vs temperature
figure('Name', 'stroke');
hold on;
plot(T_cell{5}, abs(angle_cell{5}), 'LineWidth', 1, 'Color', [  0  114  189] / 255);
plot(T, (delta(9) - l(9, :)) * 1e3 / 5.2 / pi * 180, '--', 'LineWidth', 1, 'Color', [  0  114  189] / 255);

plot(T_cell{4}, abs(angle_cell{4}), 'LineWidth', 1, 'Color', [237  177   32] / 255);
plot(T, (delta(8) - l(8, :)) * 1e3 / 5.2 / pi * 180, '--', 'LineWidth', 1, 'Color', [237  177   32] / 255);

plot(T_cell{3}, abs(angle_cell{3}), 'LineWidth', 1, 'Color', [126   47  142] / 255);
plot(T, (delta(7) - l(7, :)) * 1e3 / 5.2 / pi * 180, '--', 'LineWidth', 1, 'Color', [126   47  142] / 255);
hold off;
% legend('40 mm (Expt.)', '40 mm (Model)', '35 mm (Expt.)', '35 mm (Model)', '30 mm (Expt.)', '30 mm (Model)', 'Location', 'northwest');
xlabel('$T_{s}$\,\,($^{\circ}$C)', 'interpreter', 'latex', 'Fontsize', 8);
ylabel('$\theta$\,\,(deg)', 'interpreter', 'latex', 'Fontsize', 8);
xlim([20 60]);
ylim([0 160]);
xticks(20 : 10 : 60);
yticks(0 : 40 : 160);
grid on;
box on;
set(gca, 'FontName', 'Times New Roman','Fontsize', 8);

set(gcf, 'Units', 'Inches', 'Position', [0, 0, 1.75, 1.3], 'PaperUnits', 'Inches', 'PaperSize', [1.75, 1.3]);
saveas(gcf, '..\..\paper\v1.0.0\figures\stroke_vs_temperature.pdf');
filename = '..\..\paper\v1.0.0\figures\stroke_vs_temperature.tiff';
print(gcf, filename, '-dtiffn', '-r600');
