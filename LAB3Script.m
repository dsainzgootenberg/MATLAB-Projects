% Load Data
data = readtable('C:\Users\Dany SG\Documents\CWRU Classes\Fall 2024\EMAE 285\Labs\LAB 3\Report\MATLAB\LAB3DATA.csv');

% Extract Time and Temperatures
time = datetime(data.Time, 'InputFormat', 'MM/dd/yyyy HH:mm:ss.SSS');
T_cylinder = mean([data.Temperature_0, data.Temperature_1, data.Temperature_2], 2); % Average cylinder surface temp
T_air = data.Temperature_3;  % Ambient air temperature

% Define constants
D = 0.01905;  % Diameter of cylinder in meters
L = 0.127;    % Length of cylinder in meters
V = 12;       % Input voltage
I = 0.5;      % Input current
Pr = 0.73;    % Prandtl number
k = 0.026;    % Thermal conductivity of air (W/m-K)
epsilon = 0.1; % Emissivity of cylinder surface
sigma = 5.67e-8; % Stefan-Boltzmann constant (W/m^2-K^4)

% Define velocities in m/s
velocities = [15, 18, 21, 24];

% Initialize vectors for results
Nusselt_numbers = [];
Reynolds_numbers = [];
Q_radiation = [];
h_values = [];
errors_Nu = [];
errors_Re = [];

% Determine when the fan is running using temperature drops
fan_on_indices = find(diff(T_cylinder) < -0.1);
fan_on_periods = [fan_on_indices(1), fan_on_indices(diff(fan_on_indices) > 1)'];

% Loop over each velocity
for i = 1:length(velocities)
    % Extract data for current velocity
    if i <= length(fan_on_periods)
        start_idx = fan_on_periods(i);
        if i < length(fan_on_periods)
            end_idx = fan_on_periods(i+1) - 1;
        else
            end_idx = length(T_cylinder);
        end
    else
        break;
    end
    
    % Calculate average temperatures for this segment
    T_cyl_avg = mean(T_cylinder(start_idx:end_idx)); 
    T_air_avg = mean(T_air(start_idx:end_idx));     
    Q_conv = V * I;  % Convective heat input
    
    % Calculate heat transfer coefficient 'h'
    h = Q_conv / (pi * D * L * (T_cyl_avg - T_air_avg));
    
    % Calculate radiation heat transfer
    Q_rad = epsilon * sigma * pi * D * L * ((T_cyl_avg + 273.15)^4 - (T_air_avg + 273.15)^4);
    
    % Calculate Reynolds and Nusselt numbers
    U = velocities(i);  % Velocity for current segment
    nu = 1.5e-5;        % Kinematic viscosity of air in m^2/s
    Re = (U * D) / nu;
    Nu = (h * D) / k;
    
    % Estimate uncertainties
    err_U = 0.02 * U;
    err_Tc = 0.1 * T_cyl_avg;
    err_Ta = 0.02 * T_air_avg;
    err_Q = 0.03 * Q_conv;
    
    % Propagate uncertainties
    err_Re = Re * (err_U / U);
    err_Nu = Nu * sqrt((err_Q / Q_conv)^2 + (err_Tc / (T_cyl_avg - T_air_avg))^2);
    
    % Store results
    Reynolds_numbers = [Reynolds_numbers, Re];
    Nusselt_numbers = [Nusselt_numbers, Nu];
    Q_radiation = [Q_radiation, Q_rad];
    h_values = [h_values, h];
    errors_Re = [errors_Re, err_Re];
    errors_Nu = [errors_Nu, err_Nu];
end

% Plot Nusselt vs Reynolds Number with error bars
figure;
errorbar(Reynolds_numbers, Nusselt_numbers, errors_Nu, 'o-', 'LineWidth', 1.5);
xlabel('Reynolds Number (Re)');
ylabel('Nusselt Number (Nu)');
title('Nusselt Number vs Reynolds Number');
grid on;

% Plot Radiation vs Convection Heat Transfer
figure;
plot(velocities, Q_radiation, 'r-o', 'LineWidth', 1.5);
hold on;
plot(velocities, V * I * ones(size(velocities)), 'b--', 'LineWidth', 1.5);
xlabel('Velocity (m/s)');
ylabel('Heat Transfer Rate (W)');
legend('Radiation', 'Convection');
title('Radiation vs Convection Heat Transfer');
grid on;

% Additional plot for Temperature vs Time showing fan operation periods
figure;
plot(time, T_cylinder, 'b-', 'LineWidth', 1.5); % Cylinder surface temperature
hold on;
plot(time, T_air, 'r-', 'LineWidth', 1.5);      % Air temperature
scatter(time(fan_on_indices), T_cylinder(fan_on_indices), 50, 'r', 'filled'); % Mark fan onset points
yline(mean(T_air), '--', 'Air Temperature Average', 'LabelVerticalAlignment', 'bottom');

% Plot labels and legend
xlabel('Time');
ylabel('Temperature (°C)');
title('Cylinder Surface and Air Temperature over Time');
legend('Cylinder Temperature', 'Air Temperature', 'Fan Onset Points');
grid on;
hold off;

% Display calculated values
disp('Reynolds Numbers:'), disp(Reynolds_numbers);
disp('Nusselt Numbers:'), disp(Nusselt_numbers);
disp('Radiation Heat Transfer:'), disp(Q_radiation);
