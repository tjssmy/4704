% MATLAB script to generate animated GIF showing convergence of analytical solution for Gaussian boundary
clear all; close all; clc;

% Parameters
a = 1;              % Side length of the square
Nx = 100;           % Number of points in x for plotting
Ny = 100;           % Number of points in y for plotting
x = linspace(0, a, Nx);
y = linspace(0, a, Ny);
[X, Y] = meshgrid(x, y);

% Boundary condition at y = a
bc_top = exp( -(x - a/2).^2 / (a/5)^2 );

% Precompute coefficients for odd n only
max_modes = 100; % Maximum number of odd modes (n=1,3,...,199)
odd_ns = 1:2:(2*max_modes - 1);
B = zeros(max_modes, 1);
for ii = 1:max_modes
    n = odd_ns(ii);
    integrand = @(x_var) exp( -(x_var - a/2).^2 / (a/5)^2 ) .* sin(n * pi * x_var / a);
    int_val = integral(integrand, 0, a);
    B(ii) = (2 / a) * int_val / sinh(n * pi);
end

% Prepare for animation
fig = figure;
fig.Position(3) = fig.Position(3); % Make the figure twice as wide

% Steps for animation: add modes incrementally (1, 2, ..., max_modes)
frame_steps = 1:1:max_modes; % Number of modes per frame
num_frames = length(frame_steps);

% Create animated GIF
filename = 'convergence_analytical.gif';

for ff = 1:num_frames
    num_current_modes = frame_steps(ff);
    
    % Compute partial sum
    u_partial = zeros(size(X));
    for ii = 1:num_current_modes
        n = odd_ns(ii);
        u_partial = u_partial + B(ii) * sin(n * pi * X / a) .* sinh(n * pi * Y / a);
    end
    
    % Plot the partial solution
    % subplot(1, 2, 1);
    surf(X, Y, u_partial);
    title(sprintf('Analytical Solution (First %d Odd Modes)', num_current_modes));
    xlabel('x'); ylabel('y'); zlabel('u');
    colormap('parula'); colorbar;
    shading interp;
    hold on;
    plot3(x, a*ones(size(x)), bc_top, 'r-', 'LineWidth', 2);
    view(60, 30);
    hold off;
    
    % % Plot the boundary condition alone on the right for reference
    % subplot(1, 2, 2);
    % plot(x, bc_top, 'b-', 'LineWidth', 2);
    % title('Boundary Condition at y = a');
    % xlabel('x'); ylabel('u(x, a)');
    % ylim([0, 1.1]);
    % grid on;
    % 
    % Capture the frame
    drawnow;
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF
    if ff == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
    end
end

% Close the figure
close(fig);