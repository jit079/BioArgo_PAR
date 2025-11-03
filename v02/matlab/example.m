clc
clear

% Read example data
% 1st column: depth (m)
% 2nd column: PAR (micro E/m^2/s)
% 3rd-6th columns: Ed380, Ed443, Ed490, Ed555 (mW/cm^2/micron)
t = readtable('example/example_data.csv','VariableNamingRule', 'preserve');
df = table2array(t);

% Extract columns
d = df(:, 1); % Depth
par = df(:, 2); % Measured PAR
profile = df(:, [3, 4, 5, 6, 1]); % [Ed380, Ed443, Ed490, Ed555, depth]

% Reconstruct PAR and uncertainty
[par_est, uncertainty] = reconstruct_par(profile);

% Display estimated PAR
disp('Estimated PAR:');
disp(par_est);

% Sort data by depth for plotting
[ds, idx] = sort(d);
par_est_sorted = par_est(idx);
par_sorted = par(idx);
uncertainty_sorted = uncertainty(idx);

% Create figure with two subplots
figure('Position', [100, 100, 1000, 400]);

% Subplot 1: Measured vs Modeled PAR
subplot(1, 2, 1);
plot(ds, par_est_sorted, 'b-', 'LineWidth', 1.5, 'DisplayName', 'modeled');
hold on;
plot(ds, par_sorted, 'r-', 'LineWidth', 1.5, 'DisplayName', 'measured');
hold off;
xlabel('Depth, m', 'FontSize', 14);
ylabel('PAR, \muE/m^2/s', 'FontSize', 14);
set(gca, 'FontSize', 12, 'YScale', 'log');
legend('FontSize', 12, 'Location', 'best');

% Subplot 2: Uncertainty
subplot(1, 2, 2);
plot(ds, uncertainty_sorted, 'k-', 'LineWidth', 1.5);
xlabel('Depth, m', 'FontSize', 14);
ylabel('\DeltaPAR, \muE/m^2/s', 'FontSize', 14);
set(gca, 'FontSize', 12);

% Save figure
print('-dpng', 'example/example_measured_vs_modeled_matlab_version.png');
