% MATLAB script to solve Laplace's equation in a 2D rectangle
% using analytical (separation of variables), semi-analytical (series expansion),
% and boundary element method (BEM). Saves figures as PNG files.

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

% --- 2. Semi-Analytical: Series Expansion with Analytical Coefficients ---
% Boundary condition: u(x, b) = x*(a-x)
N = 20; % Number of series terms
u_series = zeros(ny, nx);
for n = 1:2:N % Only odd n contribute
    % Analytical Fourier coefficient for f(x) = x(a-x)
    cn = 8*(4*a^2) / (n^3 * pi^3 * sinh(n*pi*b/a));
    for i = 1:nx
        for j = 1:ny
            u_series(j, i) = u_series(j, i) + ...
                cn * sin(n*pi*x(i)/a) * sinh(n*pi*y(j)/a);
        end
    end
end

% --- 3. Semi-Analytical: Boundary Element Method (BEM) ---
% Discretize all four boundaries
N_bem = 50; % Number of points per boundary side
% Boundary points: top (y=b), bottom (y=0), left (x=0), right (x=a)
x_top = linspace(0, a, N_bem); y_top = b * ones(1, N_bem);
x_bot = linspace(0, a, N_bem); y_bot = zeros(1, N_bem);
y_left = linspace(0, b, N_bem); x_left = zeros(1, N_bem);
y_right = linspace(0, b, N_bem); x_right = a * ones(1, N_bem);

% Combine boundary points
x_b = [x_top, x_bot, x_left, y_right];
y_b = [y_top, y_bot, y_left, x_right];
N_total = 4 * N_bem;

% Boundary conditions: u on boundary
u_b = zeros(1, N_total);
u_b(1:N_bem) = sin(pi * x_top / a); % u(x, b) = sin(pi*x/a)
% u = 0 on other boundaries (x_bot, x_left, y_right)

% Normal vectors (outward): [nx, ny]
normals = zeros(N_total, 2);
normals(1:N_bem, :) = [zeros(N_bem, 1), ones(N_bem, 1)]; % Top: n = [0, 1]
normals(N_bem+1:2*N_bem, :) = [zeros(N_bem, 1), -ones(N_bem, 1)]; % Bottom: n = [0, -1]
normals(2*N_bem+1:3*N_bem, :) = [-ones(N_bem, 1), zeros(N_bem, 1)]; % Left: n = [-1, 0]
normals(3*N_bem+1:end, :) = [ones(N_bem, 1), zeros(N_bem, 1)]; % Right: n = [1, 0]

% Segment lengths
ds = [a/(N_bem-1) * ones(1, N_bem), a/(N_bem-1) * ones(1, N_bem), ...
      b/(N_bem-1) * ones(1, N_bem), b/(N_bem-1) * ones(1, N_bem)];

% Green's function and normal derivative
G = @(x, y, xp, yp) -(1/(2*pi)) * log(sqrt((x - xp).^2 + (y - yp).^2 + 1e-10));
dGdn = @(x, y, xp, yp, nx, ny) -(1/(2*pi)) * ((x - xp)*nx + (y - yp)*ny) ./ ...
    ((x - xp).^2 + (y - yp).^2 + 1e-10);

% Set up linear system for normal derivatives
A = zeros(N_total, N_total);
b = zeros(N_total, 1);
for i = 1:N_total
    for j = 1:N_total
        if i == j
            A(i, j) = 0.5; % Diagonal term for singularity
        else
            A(i, j) = -dGdn(x_b(i), y_b(i), x_b(j), y_b(j), normals(j, 1), normals(j, 2)) * ds(j);
        end
    end
    % Right-hand side
    for j = 1:N_total
        b(i) = b(i) + G(x_b(i), y_b(i), x_b(j), y_b(j)) * u_b(j) * ds(j);
    end
end

% Solve for normal derivatives
dudn = A \ b;

% Compute interior solution
u_bem = zeros(ny, nx);
for i = 1:nx
    for j = 1:ny
        sum_G = 0;
        sum_dGdn = 0;
        for k = 1:N_total
            sum_G = sum_G + G(x(i), y(j), x_b(k), y_b(k)) * dudn(k) * ds(k);
            sum_dGdn = sum_dGdn + u_b(k) * dGdn(x(i), y(j), x_b(k), y_b(k), normals(k, 1), normals(k, 2)) * ds(k);
        end
        u_bem(j, i) = sum_G - sum_dGdn;
    end
end

% --- Visualization ---
% First figure: 3D surface plots
fig1 = figure('Position', [100, 100, 1200, 400]);

% Analytical Solution
subplot(1, 3, 1);
surf(X, Y, u_analytical);
title('Analytical: Separation of Variables');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
colorbar;

% Series Expansion
subplot(1, 3, 2);
surf(X, Y, u_series);
title('Semi-Analytical: Series Expansion');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
colorbar;

% Boundary Element Method
subplot(1, 3, 3);
surf(X, Y, u_bem);
title('Semi-Analytical: Boundary Element Method');
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
colorbar;

sgtitle('Solutions to Laplace''s Equation');

% Save first figure as PNG
saveas(fig1, 'laplace_solutions.png');

% --- Validation for BEM ---
% Check boundary condition u(x, b) = sin(pi*x/a) at y = b
% x_check = linspace(0, a, 100);
% u_check = zeros(size(x_check));
% for k = 1:length(x_check)
%     sum_G = 0;
%     sum_dGdn = 0;
%     for m = 1:N_total
%         sum_G = sum_G + G(x_check(k), b, x_b(m), y_b(m)) * dudn(m) * ds(m);
%         sum_dGdn = sum_dGdn + u_b(m) * dGdn(x_check(k), b, x_b(m), y_b(m), normals(m, 1), normals(m, 2)) * ds(m);
%     end
%     u_check(k) = sum_G - sum_dGdn;
% end
% fig2 = figure;
% plot(x_check, sin(pi*xhammer_check/a), 'b-', 'DisplayName', 'Exact: sin(\pi x/a)');
% hold on;
% plot(x_check, u_check, 'r--', 'DisplayName', 'BEM at y=b');
% xlabel('x'); ylabel('u(x, b)'); title('BEM Solution at y = b');
% legend; grid on;
% 
% % Save second figure as PNG
% saveas(fig2, 'bem_validation.png');