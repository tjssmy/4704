set(0,'defaultaxesfontname','Times New Roman');
set(0,'DefaultLineLineWidth', 2);
% set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultaxesfontsize',20);



% Constants
m = 9.11e-31;       % Electron mass (kg)
hbar = 1.0545718e-34; % Reduced Planck constant (J·s)
eV_to_J = 1.60218e-19; % eV to Joules

% Fixed energy
E = 5 * eV_to_J;    % Particle energy, 5 eV in Joules

% Variable ranges
V0_eV = linspace(6, 15, 50); % Barrier height in eV (6 to 15 eV)
V0 = V0_eV * eV_to_J;        % Convert to Joules
F = logspace(-12, -9, 50);   % Force in N (10^-12 to 10^-9 N)

% Fixed parameters for 1D plots
F_fixed = 1e-10;    % Fixed force for T vs V0 plot (N)
V0_fixed = 10 * eV_to_J; % Fixed barrier height for T vs F plot (J)

% Calculate T vs V0
T_V0 = exp(-4 * sqrt(2 * m) * (V0 - E).^(3/2) / (3 * hbar * F_fixed));

% Calculate T vs F
T_F = exp(-4 * sqrt(2 * m) * (V0_fixed - E)^(3/2) ./ (3 * hbar * F));

% Meshgrid for surface plot
[V0_mesh, F_mesh] = meshgrid(V0, F);
T_surface = exp(-4 * sqrt(2 * m) * (V0_mesh - E).^(3/2) ./ (3 * hbar * F_mesh));

% Plot 1: T vs V0
figure(1);
semilogy(V0_eV, T_V0, 'b-', 'LineWidth', 2);
xlabel('Barrier Height V_0 (eV)');
ylabel('Tunneling Probability T');
title('T vs V_0 (E = 5 eV, F = 10^{-10} N)');
grid on;
ylim([1e-10  1e-5]);
saveas(gcf, 'T_vs_V0.png');

% Plot 2: T vs F
figure(2);
loglog(F, T_F, 'r-', 'LineWidth', 2);
xlabel('Force F (N)');
ylabel('Tunneling Probability T');
title('T vs F (E = 5 eV, V_0 = 10 eV)');
grid on;
ylim([1e-10 1e-5]);
saveas(gcf, 'T_vs_F.png');

% Plot 3: Surface plot of T vs V0 and F
figure(3);
surf(V0_eV, F, log10(T_surface), 'EdgeColor', 'none');
shading interp
set(gca, 'YScale', 'log');
xlabel('V_0 (eV)');
ylabel('F (N)');
zlabel('T');
title('T vs V_0 and F (E = 5 eV)');
colorbar;
zlim([-50 0]);
clim([-50 0]);
saveas(gcf, 'T_vs_V0_and_F.png');