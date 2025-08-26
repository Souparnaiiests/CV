clc; clear; close all;

%% INPUT PARAMETERS
H = 10;               % Thickness of soil layer (m)
Nz = 100;             % Number of spatial divisions
dz = H / Nz;          % Grid spacing
cv = 8.6881e-5;       % Coefficient of consolidation (m^2/s)

% Stability condition for implicit scheme
lambda_max = 0.5;
dt_max = lambda_max * dz^2 / cv;
dt = dt_max * 0.9;     % Safe time step for stability
T = 3600 * 24 * 10;    % Total simulation time (s)
Nt = round(T / dt);    % Number of time steps
lambda = cv * dt / dz^2;

% Semi-permeable boundary resistances (dimensionless)
Rt = 10;    % Top boundary resistance
Rb = 10;    % Bottom boundary resistance

sigma0 = 3.47;                % Initial effective stress (kPa)
sigma_dash = 16 * sigma0;     % Final effective stress (kPa)

%% GRID & INITIAL CONDITIONS
z = linspace(0, H, Nz+1);     % Depth grid
t = linspace(0, T, Nt+1);     % Time grid
w = zeros(Nz+1, Nt+1);        % Initialize transformed variable w
w(:,1) = log10(sigma_dash / sigma0);

%% ASSEMBLE MATRIX (IMPLICIT SCHEME WITH SEMI-PERMEABLE BOUNDARIES)
main_diag = (1 + 2*lambda) * ones(Nz+1, 1);
off_diag = -lambda * ones(Nz, 1);
A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

% Modify top boundary at z = 0
A(1,1) = 1 + lambda + lambda * dz * Rt;
A(1,2) = -lambda;

% Modify bottom boundary at z = H
A(end,end) = 1 + lambda + lambda * dz * Rb;
A(end,end-1) = -lambda;

%% TIME MARCHING SOLVER
for n = 1:Nt
    b = w(:,n);
    w_new = A \ b;
    w(:,n+1) = w_new;
end

%% POST-PROCESSING
u = sigma_dash * (1 - 10.^(-w));        % Recover pore pressure
u_max = max(u);                         % Maximum excess pore pressure over depth
u_max_percent = (u_max / u_max(1)) * 100; % Normalize as percentage
Tv = cv * t / (H)^2;                   % Time factor (T_v)

%% PLOT RESULTS
figure;
semilogx(Tv, u_max_percent, 'b-', 'LineWidth', 2); hold on;
xlabel('Time factor T_v');
ylabel('Maximum Excess Pore Pressure (\% of Initial)');
title('Numerical Results: Semi-Permeable Boundaries');
legend('Semi-Permeable (This Study)', 'Location', 'northeast');
ylim([0 110]);
grid on;