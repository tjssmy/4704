% Analytical solution to 1D time-dependent heat equation using Fourier series
% Boundary conditions: T(0,t)=0, T(L,t)=0
% Initial condition: T(x,0) = exp( -(x - L/3)^2 / (2*(L/8)^2) )

clear; close all; clc;

% Parameters
L = 1;              % Length of the domain
alpha = 0.1;        % Thermal diffusivity
N = 100;            % Number of Fourier terms (increase for better accuracy)
Nx = 200;           % Number of spatial points for plotting
Nt = 100;           % Number of time steps for animation
t_max = 10;          % Maximum time for simulation

% Spatial grid
x = linspace(0, L, Nx);

% Time grid
t = linspace(0, t_max, Nt);

% Precompute eigenvalues and eigenfunctions
lambda_n = @(n) (n * pi / L)^2;
X_n = @(n, x) sin(n * pi * x / L);

% Initial condition function (Gaussian)
sigma = L / 8;
f = @(x) exp( -(x - L/3).^2 / (2 * sigma^2) );

% Compute Fourier coefficients c_n
c_n = zeros(1, N);
for n = 1:N
    integrand = @(x) f(x) .* X_n(n, x);
    c_n(n) = (2 / L) * integral(integrand, 0, L);
end

% Prepare figure for animation
figure;
hold on;
plot(x, f(x), 'r--', 'LineWidth', 2); % Initial temperature as dotted red line
h = plot(x, f(x), 'b-', 'LineWidth', 2); % Evolving temperature as solid blue line
hold off;
xlabel('x');
ylabel('T(x,t)');
title('Time Evolution of Temperature');
legend('Initial T(x,0)', 'T(x,t)');
ylim([0, max(f(x)) * 1.1]);
grid on;

% GIF filename
filename = 'heat_evolution.gif';

% Capture initial frame
frame = getframe(gcf);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);

% Animate over time
for k = 2:Nt
    T = zeros(1, Nx);
    for n = 1:N
        T = T + c_n(n) * X_n(n, x) * exp(-alpha * lambda_n(n) * t(k));
    end
    set(h, 'YData', T);
    drawnow;
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
end