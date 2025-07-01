%% PETM Temperature Reconstruction - Monte Carlo Analysis
% Site U1580 Acarinina coalingensis Mg/Ca paleothermometry
% Author: Samina Anee, Utah State University
% Date: 2024

clear all; close all; clc;

%% Input Data - Acarinina coalingensis Mg/Ca values (mmol/mol)
% Pre-PETM Acarinina coalingensis Mg/Ca values
MgCa_prePETM = [2.419, 2.702, 2.075, 2.884, 2.801, 2.592]; % Pre-PETM A. coalingensis
% PETM Acarinina coalingensis Mg/Ca values  
MgCa_PETM = [3.561, 3.193, 3.626, 4.08]; % PETM A. coalingensis

%% Analysis Parameters
num_samples = 1000;          % Number of Monte Carlo iterations
salinity_error = 0.5;        % PSU uncertainty
pH_error = 0.1;              % pH uncertainty  
analytical_error = 0.03;     % 3% analytical uncertainty on Mg/Ca

% Environmental Parameters
ph_prePETM = 8.05;           % Pre-PETM pH
pH_PETM = 7.8;               % PETM pH (ocean acidification)
salinity = 35;               % Assumed salinity (PSU)
MgCa_sw = 2.0;               % Paleocene-Eocene seawater Mg/Ca ratio (mol/mol)

% Calibration Constants (Anand et al., 2003)
H = 0.41;                    % Power law coefficient (Hasiuk & Lohmann, 2010)
B_modern = 0.38;             % Modern calibration constant for A. coalingensis
A = 0.09;                    % Temperature sensitivity coefficient

%% Prepare data arrays
MgCa_mean = [mean(MgCa_prePETM), mean(MgCa_PETM)];
MgCa_std = [std(MgCa_prePETM), std(MgCa_PETM)];

% Initialize output arrays
T_est = zeros(2, num_samples); % 1=pre-PETM, 2=PETM
MgCa_corrected = zeros(2, num_samples);

%% Monte Carlo Simulation
fprintf('Running Monte Carlo simulation with %d iterations...\n', num_samples);

for i = 1:2  % 1 = pre-PETM, 2 = PETM
    for j = 1:num_samples
        
        % Add uncertainty to raw Mg/Ca values
        MgCa_raw = normrnd(MgCa_mean(i), sqrt(MgCa_std(i)^2 + (MgCa_mean(i)*analytical_error)^2));
        
        % Apply salinity corrections with uncertainty
        sal_local = normrnd(salinity, salinity_error);
        sal_corrected = (1-(sal_local - 35)*0.042) * MgCa_raw;
        
        % Apply pH corrections based on time period
        if i == 1  % pre-PETM
            pH_local = normrnd(ph_prePETM, pH_error);
        else       % PETM
            pH_local = normrnd(pH_PETM, pH_error);
        end
        
        % Calculate final corrected Mg/Ca values
        MgCa_corrected(i,j) = (1-(8.05-pH_local)*0.7) * sal_corrected;
        
        % Apply seawater Mg/Ca correction
        B_corrected = (MgCa_sw^H / 5.2^H) * B_modern;
        
        % Temperature calculation
        T_est(i,j) = log(MgCa_corrected(i,j) / B_corrected) / A;
    end
end

%% Calculate Statistics
pre_PETM_mean = mean(T_est(1,:));
peak_PETM_mean = mean(T_est(2,:));
temp_change = peak_PETM_mean - pre_PETM_mean;

pre_PETM_std = std(T_est(1,:));
peak_PETM_std = std(T_est(2,:));
temp_change_uncertainty = sqrt(pre_PETM_std^2 + peak_PETM_std^2);

% Calculate percentiles for confidence intervals
pre_PETM_ci = prctile(T_est(1,:), [2.5, 97.5]);
PETM_ci = prctile(T_est(2,:), [2.5, 97.5]);

%% Display Results
fprintf('\n=== Acarinina coalingensis PETM Temperature Reconstruction Results ===\n');
fprintf('Pre-PETM Temperature (A. coalingensis): %.1f +/- %.1f°C (95%% CI: %.1f-%.1f°C)\n', ...
    pre_PETM_mean, pre_PETM_std, pre_PETM_ci(1), pre_PETM_ci(2));
fprintf('Peak PETM Temperature (A. coalingensis): %.1f +/- %.1f°C (95%% CI: %.1f-%.1f°C)\n', ...
    peak_PETM_mean, peak_PETM_std, PETM_ci(1), PETM_ci(2));
fprintf('Temperature Change (A. coalingensis): %.1f +/- %.1f°C\n', temp_change, temp_change_uncertainty);
fprintf('Raw Mg/Ca change (A. coalingensis): %.2f mmol/mol\n', MgCa_mean(2) - MgCa_mean(1));

%% Create Publication-Quality Figure with Panel Labels
figure('Position', [100 100 800 600], 'Color', 'white');

% Panel (a) - Main histogram plot
subplot(2,2,[1,2]);
h1 = histogram(T_est(1,:), 'BinWidth', 0.5, 'FaceColor', [0.8 0.4 0.2], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'probability');
hold on;
h2 = histogram(T_est(2,:), 'BinWidth', 0.5, 'FaceColor', [0.2 0.6 0.8], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.7, 'Normalization', 'probability');

% Add mean lines
line([pre_PETM_mean pre_PETM_mean], ylim, 'Color', [0.8 0.4 0.2], 'LineWidth', 2, 'LineStyle', '--');
line([peak_PETM_mean peak_PETM_mean], ylim, 'Color', [0.2 0.6 0.8], 'LineWidth', 2, 'LineStyle', '--');

xlabel('Temperature (°C)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Probability Density', 'FontSize', 12, 'FontWeight', 'bold');
title('Acarinina coalingensis Mg/Ca Derived Temperature Distributions', 'FontSize', 14, 'FontWeight', 'bold');
legend([h1, h2], {'Pre-PETM', 'Peak PETM'}, 'Location', 'northwest', 'FontSize', 11);
grid on; grid minor;
xlim([20 35]);

% Add panel label (a)
text(21, 0.14, '(a)', 'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Panel (b) - Box plot comparison
subplot(2,2,3);
boxplot([T_est(1,:)', T_est(2,:)'], {'Pre-PETM', 'PETM'}, 'Colors', [0.8 0.4 0.2; 0.2 0.6 0.8]);
ylabel('Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('Temperature Comparison', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Add panel label (b)
text(0.6, 32, '(b)', 'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'black');

% Panel (c) - Mg/Ca vs Temperature relationship
subplot(2,2,4);
scatter(mean(MgCa_corrected(1,:)), pre_PETM_mean, 100, [0.8 0.4 0.2], 'filled', 'MarkerEdgeColor', 'k');
hold on;
scatter(mean(MgCa_corrected(2,:)), peak_PETM_mean, 100, [0.2 0.6 0.8], 'filled', 'MarkerEdgeColor', 'k');

% Add error bars
errorbar(mean(MgCa_corrected(1,:)), pre_PETM_mean, pre_PETM_std, pre_PETM_std, ...
    std(MgCa_corrected(1,:)), std(MgCa_corrected(1,:)), 'Color', [0.8 0.4 0.2], 'LineWidth', 2);
errorbar(mean(MgCa_corrected(2,:)), peak_PETM_mean, peak_PETM_std, peak_PETM_std, ...
    std(MgCa_corrected(2,:)), std(MgCa_corrected(2,:)), 'Color', [0.2 0.6 0.8], 'LineWidth', 2);

xlabel('Corrected Mg/Ca (mmol/mol)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Temperature (°C)', 'FontSize', 11, 'FontWeight', 'bold');
title('Mg/Ca-Temperature Relationship', 'FontSize', 12, 'FontWeight', 'bold');
legend({'Pre-PETM', 'PETM'}, 'Location', 'southeast', 'FontSize', 10);
grid on;

% Add panel label (c)
text(0.05, 0.95, '(c)', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', 'white', 'EdgeColor', 'black');
% Overall figure title
sgtitle('Acarinina coalingensis PETM Temperature Reconstruction - Site U1580', 'FontSize', 16, 'FontWeight', 'bold');

%% Save Results
% Create results table
Results = table();
Results.Period = {'Pre-PETM'; 'PETM'};
Results.Mean_Temp_C = [pre_PETM_mean; peak_PETM_mean];
Results.Std_Temp_C = [pre_PETM_std; peak_PETM_std];
Results.CI_Lower_C = [pre_PETM_ci(1); PETM_ci(1)];
Results.CI_Upper_C = [pre_PETM_ci(2); PETM_ci(2)];
Results.Mean_MgCa_mmol_mol = [mean(MgCa_corrected(1,:)); mean(MgCa_corrected(2,:))];

fprintf('\n=== Results Table ===\n');
disp(Results);

% Save figure and data
saveas(gcf, 'PETM_Acarinina_coalingensis_Monte_Carlo.png', 'png');
save('PETM_Acarinina_coalingensis_Results.mat', 'T_est', 'MgCa_corrected', 'Results');

%% Additional Analysis: Probability of Warming
warming_probability = sum(T_est(2,:) > T_est(1,:)) / num_samples * 100;
fprintf('\nProbability of PETM warming (A. coalingensis): %.1f%%\n', warming_probability);

% Temperature difference distribution
temp_diff = T_est(2,:) - T_est(1,:);
fprintf('Temperature difference (A. coalingensis): %.1f +/- %.1f°C\n', mean(temp_diff), std(temp_diff));

fprintf('\nAcarinina coalingensis analysis complete! Figure and data saved.\n');