% MATLAB script to solve 1D Laplace's equation with finite differences
% Boundary conditions: du/dx(0) = 2, u(L) = 0
% Corrected analytical solution
% Saves figures as PNG files

clear all;
close all;

% Parameters
L = 1; % Domain length
N = 50; % Number of intervals
dx = L / N; % Grid spacing
x = 0:dx:L; % Grid points

% Initialize coefficient matrix A and RHS vector b
A = zeros(N+1, N+1);
b = zeros(N+1, 1);

% Interior points: u_{i+1} - 2u_i + u_{i-1} = 0
for i = 2:N
    A(i, i-1) = 1;
    A(i, i) = -2;
    A(i, i+1) = 1;
    b(i) = 0;
end

% Neumann BC at x=0: (u_1 - u_0)/dx = 2
A(1, 1) = -1;
A(1, 2) = 1;
b(1) = 2 * dx;

% Dirichlet BC at x=L: u_N = 0
A(N+1, N+1) = 1;
b(N+1) = 0;

% Solve the system
u = A \ b;

% Corrected analytical solution: u(x) = 2x - 2L
u_analytical = 2 * x - 2 * L;

% Plot 1: Geometry and Discretization
figure(1);
plot(x, zeros(size(x)), 'b-', 'LineWidth', 2);
hold on;
plot(x, zeros(size(x)), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
plot(0, 0, 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(L, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
text(-0.05, 0.1, 'du/dx = 2', 'FontSize', 12, 'Color', 'g');
text(L-0.05, 0.1, 'u = 0', 'FontSize', 12, 'Color', 'k');
xlabel('x');
ylabel('y');
title('Geometry and Discretization');
legend('Domain [0, L]', 'Grid points', 'Neumann BC', 'Dirichlet BC', 'Location', 'NorthWest');
grid on;
axis([-0.1 L+0.1 -0.5 0.5]);
hold off;
% Save figure as PNG
print('-dpng', 'geometry_discretization.png');

% Plot 2: Solution
figure(2);
plot(x, u, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Numerical');
hold on;
plot(x, u_analytical, 'r--', 'LineWidth', 2, 'DisplayName', 'Analytical: u(x) = 2x - 2L');
xlabel('x');
ylabel('u(x)');
title('Solution of 1D Laplace''s Equation');
legend('show');
grid on;
hold off;
% Save figure as PNG
print('-dpng', 'solution.png');