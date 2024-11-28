% MATLAB Script for Vortex Panel Method - Lab 4 Airfoil Experiment
% Calculates lift and moment coefficients and ensures plots for specified AOAs.

% Clear workspace and command window
clear;
clc;
close all;

% Load experimental data
data = readmatrix('C:\\Users\\Dany SG\\Documents\\CWRU Classes\\Fall 2024\\EMAE 285\\Labs\\LAB 4\\Report\\MATLAB\\LAB_4_dataCSV.csv');
AOA = data(:, 4);          % Angle of attack (degrees)
TunnelVelocity = data(:, 2); % Tunnel velocity (m/s)
P_US = data(:, 21:28);     % Upper surface pressure readings
P_LS = data(:, 29:36);     % Lower surface pressure readings

% Constants
rho = 1.225; % Air density (kg/m^3)
P_inf = 0;   % Free-stream pressure (Pa), assumed zero-gauge
chord = 1;   % Normalized chord length

% Generate chord positions for panels
numPoints = size(P_US, 2); % Number of pressure tap points
x = linspace(0, 1, numPoints); % Chord-wise positions

% Calculate pressure coefficients (Cp)
Cp_US = (P_US - P_inf) ./ (0.5 * rho * TunnelVelocity.^2); % Upper surface
Cp_LS = (P_LS - P_inf) ./ (0.5 * rho * TunnelVelocity.^2); % Lower surface

% Preallocate arrays for lift and moment coefficients
LiftCoefficient = zeros(size(AOA));
MomentCoefficient = zeros(size(AOA));

% Specify angles of attack for which to plot Cp
selectedAOAs = [-10, 0, 10];

% Loop through each angle of attack
for k = 1:length(AOA)
    % Calculate Cp difference
    Cp_diff = Cp_LS(k, :) - Cp_US(k, :);

    % Lift coefficient (C_L)
    LiftCoefficient(k) = trapz(x, Cp_diff);

    % Moment coefficient (C_M)
    % Assume moment arm about the quarter-chord point
    moment_arm = x - 0.25; % Distance from quarter chord
    MomentCoefficient(k) = trapz(x, moment_arm .* Cp_diff);

    % Find and plot for the closest AOAs (-10°, 0°, 10°)
    for i = 1:length(selectedAOAs)
        if abs(AOA(k) - selectedAOAs(i)) < 0.1 % Allow a small tolerance for matching
            figure;
            plot(x, Cp_US(k, :), '-o', 'DisplayName', 'Upper Surface');
            hold on;
            plot(x, Cp_LS(k, :), '-s', 'DisplayName', 'Lower Surface');
            xlabel('Chord Position');
            ylabel('Pressure Coefficient (Cp)');
            title(sprintf('Pressure Coefficient Distribution (AOA = %.1f°)', AOA(k)));
            legend;
            grid on;
        end
    end
end

% Plot lift coefficient vs. angle of attack
figure;
plot(AOA, LiftCoefficient, '-x');
xlabel('Angle of Attack (degrees)');
ylabel('Lift Coefficient (C_L)');
title('Lift Coefficient vs. Angle of Attack');
grid on;

% Plot moment coefficient vs. angle of attack
figure;
plot(AOA, MomentCoefficient, '-o');
xlabel('Angle of Attack (degrees)');
ylabel('Moment Coefficient (C_M)');
title('Moment Coefficient vs. Angle of Attack');
grid on;

% Output results
fprintf('Lift Coefficient Values:\n');
disp(LiftCoefficient);
fprintf('Moment Coefficient Values:\n');
disp(MomentCoefficient);
