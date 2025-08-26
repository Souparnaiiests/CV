clc
clear
close all

% Parameters
H = 10.0;                % Thickness of soil layer (m)
Nz = 100;                % Number of spatial divisions
dz = H / Nz;
cv = 1e-4;              % Coefficient of consolidation (m^2/s)
dt = 10;                 % Time step (s)
T = 3600*50;               % Total time (s)
Nt = round(T / dt);
lambda = cv * dt / dz^2;

% Grid and time
z = linspace(0, H, Nz+1);
t = linspace(0, T, Nt+1);

% Initial condition
u = zeros(Nz+1, Nt+1);   % u(z,t)
u(:,1) = 100;            % Initial excess pore pressure (kPa)

% Boundary conditions: double drainage (u=0 at top and bottom)
u(1,:) = 0;
u(end,:) = 0;

% Coefficient matrix A for implicit scheme (interior nodes only)
main_diag = (1 + 2*lambda) * ones(Nz-1, 1);
off_diag = -lambda * ones(Nz-2, 1);
A = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

% Time stepping loop
for n = 1:Nt
    % Right-hand side is u(:,n) at interior nodes
    b = u(2:Nz, n);

    % Solve A * u_new = b
    u_new = A \ b;

    % Update u (keep boundaries at 0)
    u(2:Nz, n+1) = u_new;
end

% Compute Uavg and Tv

Tv = cv * t / (H/2)^2;
area(1)=100*H;

for n = 1:Nt+1
    area(n)=polyarea(u(:,n),z');
    Uavg(n) = 1 - (area(n)/area(1));
end

% Theoretical solution for comparison
Uavg_theory=0.01:0.01:0.99; % take 

for ii=1:1:length(Uavg_theory)
    if Uavg_theory(ii)<=0.6
        Tv_theory(ii)=(pi/4).*(Uavg_theory(ii)^2);
    else
        % tv2(ii)=1.7781-0.933*log10(100-(Uavg_theory(ii)*100));
        Tv_theory(ii)=-0.085-0.933*log10(1-(Uavg_theory(ii)));
    end
end





% Plot Uavg vs Tv
figure;
semilogx(Tv, Uavg, 'b-', 'LineWidth', 2);
hold on;
semilogx(Tv_theory, Uavg_theory, 'r--', 'LineWidth', 2);
hold on;
set(gca,'Ydir','reverse');
xlabel('Time factor T_v');
ylabel('Average Degree of Consolidation U_{avg}');
title('U_{avg} vs T_v: Implicit Scheme vs Theoretical');
legend('Numerical (Implicit)', 'Theoretical');
grid on;

% Plot final pressure profile
figure;
plot(u(:,end), z, 'k-', 'LineWidth', 2);
set(gca, 'YDir', 'reverse');
xlabel('Excess pore pressure (kPa)');
ylabel('Depth (m)');
title('Pore Pressure Distribution at Final Time (Implicit Scheme)');
grid on;
