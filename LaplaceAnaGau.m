% MATLAB script to solve 2D Laplace equation on a square domain with Gaussian boundary condition
clear all; close all; clc;

% Parameters
a = 1;              % Side length of the square
N = 50;             % Number of grid points in each direction (N+2 includes boundaries)
dx = a / (N+1);     % Grid spacing
x = 0:dx:a;         % Grid points in x
y = 0:dx:a;         % Grid points in y
[X, Y] = meshgrid(x, y);

% Initialize solution matrix (including boundaries)
u = zeros(N+2, N+2);

% Apply boundary conditions
% u(0, y) = 0 (left), u(a, y) = 0 (right), u(x, 0) = 0 (bottom) are zero by initialization
% u(x, a) = exp( -(x - a/2)^2 / (a/5)^2 ) (top)
u(N+2, :) = exp( -(x - a/2).^2 / (a/5)^2 );

% Set up the linear system for interior points
N_interior = N * N; % Number of interior points
A = sparse(N_interior, N_interior); % Coefficient matrix
b = zeros(N_interior, 1); % Right-hand side vector

% Map 2D grid to 1D vector: (i,j) -> k = i + (j-1)*N
for j = 1:N
    for i = 1:N
        k = i + (j-1)*N; % Index in the 1D vector
        % Finite difference for Laplace: (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j))/dx^2 = 0
        A(k, k) = -4 / dx^2; % Central point
        % Neighboring points
        if i > 1
            A(k, (i-1) + (j-1)*N) = 1 / dx^2; % u(i-1,j)
        else
            b(k) = b(k) - u(1, j+1) / dx^2; % Boundary at x=0
        end
        if i < N
            A(k, (i+1) + (j-1)*N) = 1 / dx^2; % u(i+1,j)
        else
            b(k) = b(k) - u(N+2, j+1) / dx^2; % Boundary at x=a
        end
        if j > 1
            A(k, i + (j-2)*N) = 1 / dx^2; % u(i,j-1)
        else
            b(k) = b(k) - u(i+1, 1) / dx^2; % Boundary at y=0
        end
        if j < N
            A(k, i + j*N) = 1 / dx^2; % u(i,j+1)
        else
            b(k) = b(k) - u(i+1, N+2) / dx^2; % Boundary at y=a
        end
    end
end

% Solve the linear system
u_vec = A \ b;

% Reshape solution to 2D grid (interior points only)
u_interior = reshape(u_vec, [N, N]);
u(2:N+1, 2:N+1) = u_interior;

% Compute analytical solution using Fourier series (only odd n due to symmetry)
N_modes = 100; % Number of odd modes to include
odd_ns = 1:2:(2*N_modes - 1);
B = zeros(N_modes, 1);
for ii = 1:N_modes
    n = odd_ns(ii);
    int_val = integral(@(x_var) exp( -(x_var - a/2).^2 / (a/5)^2 ) .* sin(n * pi * x_var / a), 0, a);
    B(ii) = (2 / a) * int_val / sinh(n * pi);
end

u_analytical = zeros(size(X));
for ii = 1:N_modes
    n = odd_ns(ii);
    u_analytical = u_analytical + B(ii) * sin(n * pi * X / a) .* sinh(n * pi * Y / a);
end

% Plot numerical solution
fig = figure;
fig.Position(3) = 2 * fig.Position(3); % Make the figure twice as wide
subplot(1, 2, 1);
surf(X, Y, u);
title('Numerical Solution');
xlabel('x'); ylabel('y'); zlabel('u');
colormap('parula'); colorbar;
shading interp;
hold on;
plot3(x, a*ones(size(x)), exp( -(x - a/2).^2 / (a/5)^2 ), 'r-', 'LineWidth', 2);
view(60, 30);

% Plot analytical solution
subplot(1, 2, 2);
surf(X, Y, u_analytical);
title('Analytical Solution (Series)');
xlabel('x'); ylabel('y'); zlabel('u');
colormap('parula'); colorbar;
shading interp;
hold on;
plot3(x, a*ones(size(x)), exp( -(x - a/2).^2 / (a/5)^2 ), 'r-', 'LineWidth', 2);
view(60, 30);

% Save the figure as PNG
print('-dpng', 'laplace_solutions_gaussian.png');

% Compute and display maximum absolute error
error = max(max(abs(u - u_analytical)));
fprintf('Maximum absolute error: %e\n', error);