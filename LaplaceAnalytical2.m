% MATLAB script to solve Laplace's equation in a 2D rectangle
% using analytical (separation of variables), semi-analytical (series expansion),
% and boundary element method (BEM).

clear all; close all;

% Domain parameters
a = 1; % width of rectangle
b = 1; % height of rectangle
nx = 50; % grid points in x
ny = 50; % grid points in y
x = linspace(0, a, nx);
y = linspace(0, b, ny);
[X, Y] = meshgrid(x, y);

% --- 1. Analytical Solution: Separation of Variables ---
% Boundary condition: u(x, b) = sin(pi*x/a)
u_analytical = zeros(ny, nx);
n = 1; % Only n=1 term for f(x) = sin(pi*x/a)
for i = 1:nx
    for j = 1:ny
        u_analytical(j, i) = (sin(pi*x(i)/a) * sinh(pi*y(j)/a)) / sinh(pi*b/a);
    end
end

% --- 2. Semi-Analytical: Series Expansion with Numerical Coefficients ---
% Boundary condition: u(x, b) = x*(a-x)
N = 200; % Number of series terms
u_series = zeros(ny, nx);
for n = 1:N
    % Compute Fourier coefficient numerically
    f = @(x) x.*(a-x).*sin(n*pi*x/a);
    % f = @(x) x.*(a-x);
    cn = (2/a) * integral(f, 0, a) / sinh(n*pi*b/a);
    for i = 1:nx
        for j = 1:ny
            u_series(j, i) = u_series(j, i) + ...
                cn * sin(n*pi*x(i)/a) * sinh(n*pi*y(j)/a);
        end
    end
end

% % --- 3. Semi-Analytical: Boundary Element Method (BEM) ---
% % Discretize top boundary (y = b) for simplicity
% N_bem = 50; % Number of boundary points
% x_b = linspace(0, a, N_bem);
% ds = a / (N_bem - 1); % Boundary segment length
% u_bem = zeros(ny, nx);
% G = @(r, rp) -(1/(2*pi)) * log(abs(r - rp)); % 2D Green's function
% dGdn = @(r, rp, n) -(1/(2*pi)) * (r - rp)./(abs(r - rp).^2 + 1e-10); % Normal derivative
% 
% % Solve for potential at interior points
% for i = 1:nx
%     for j = 1:ny
%         r = [x(i); y(j)];
%         sum_G = 0;
%         sum_dGdn = 0;
%         for k = 1:N_bem
%             rp = [x_b(k); b]; % Point on top boundary
%             n = [0; -1]; % Outward normal at y = b
%             % Boundary condition u(x, b) = sin(pi*x/a)
%             u_b = sin(pi*x_b(k)/a);
%             % Approximate normal derivative (simplified assumption)
%             dudn = 0; % Assume zero for Dirichlet problem simplification
%             sum_G = sum_G + G(r, rp) * dudn * ds;
%             sum_dGdn = sum_dGdn + u_b * dGdn(r, rp, n) * ds;
%         end
%         u_bem(j, i) = sum_G - sum_dGdn;
%     end
% end

% --- Visualization ---
figure('Position', [100, 100, 1200, 400]);

% Analytical Solution
subplot(1, 2, 1);
surf(X, Y, u_analytical);
title('Analytical: Separation of Variables');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
colorbar;

% Series Expansion
subplot(1, 2, 2);
surf(X, Y, u_series);
title('Semi-Analytical: Series Expansion');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
colorbar;

% % Boundary Element Method
% subplot(1, 3, 3);
% surf(X, Y, u_bem);
% title('Semi-Analytical: Boundary Element Method');
% xlabel('x'); ylabel('y'); zlabel('u(x,y)');
% colorbar;

sgtitle('Solutions to Laplace''s Equation');