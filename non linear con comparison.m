clc; clear; close all;

%% =================== INPUT PARAMETERS (Validation Case) ===================
H = 10;                    % Thickness of soil layer (m)
Nz = 100;                  % Number of spatial divisions
dz = H / Nz;
cv = 8.6881e-5;            % Coefficient of consolidation (m^2/s)

% Time-step parameters
lambda_max = 0.5;
dt_max = lambda_max * dz^2 / cv;
dt = dt_max * 0.9;         % Safe time-step
T = 3600 * 24 * 10;        % Total time (10 days)
Nt = round(T / dt);        % Number of time-steps
lambda = cv * dt / dz^2;

%% Effective Stress Parameters
sigma0 = 3.47;
sigma_dash = 16 * sigma0;

%% =================== GRID & INITIALIZATION ===================
z = linspace(0, H, Nz+1);      % Spatial grid (depth)
t = linspace(0, T, Nt+1);      % Time vector

% Initial log-transformed excess pressure
w = zeros(Nz+1, Nt+1);
w(:,1) = log10(sigma_dash / sigma0);   % Initial excess pressure (log10)

% Boundary conditions
w(1,:) = 0;                   % Top drainage boundary (u=0)
w(end-1,:) = w(end,:);        % Bottom impermeable boundary (∂u/∂z = 0)

%% =================== IMPLICIT TIME-STEPPING SOLVER ===================
main_diag = (1 + 2*lambda) * ones(Nz, 1);
off_diag = -lambda * ones(Nz-1, 1);
A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
A(Nz, Nz-1) = 2 * A(Nz, Nz-1);  % Adjust for PTIB

for n = 1:Nt
    b = w(2:Nz+1, n);               % Interior nodes
    w_new = A \ b;                 % Solve linear system
    w(2:Nz+1, n+1) = w_new;        % Update solution
end

%% =================== PORE PRESSURE & MAXIMUM % ===================
u = sigma_dash * (1 - 10.^(-w));        % Excess pore pressure (kPa)
u_max = max(u);                         % Max u at each time step
u_max_percent = (u_max / u_max(1)) * 100;

Tv = cv * t / (H^2);                    % Time factor

%% =================== EXPERIMENTAL DATA FROM PAPER ===================
Tv_exp = [0.02 0.04 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4];
u_percent_exp = [100 100 99.5 99 94.1 86.8 78.1 68.5 58.9 49.7 41.3 34.0 27.6 2.67 0.23 0.02];

%% =================== PLOT: COMPARISON ===================
figure;
semilogx(Tv, u_max_percent, 'm-', 'LineWidth', 2); hold on;
semilogx(Tv_exp, u_percent_exp, 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
xlabel('Time factor T_v');
ylabel('Maximum Excess Pore Pressure (% of Initial)');
title('Comparison of Numerical Model vs Experimental Data (σ''_f / σ''_0 = 16)');
legend('Numerical (This Study)', 'Experimental (Davis & Raymond, 1965)', 'Location', 'northeast');
ylim([0 110]); grid on;