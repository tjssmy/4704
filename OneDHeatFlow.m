% 1D Heat Equation Solution with Gaussian Source (Fourier Series)
clear; close all; clc;

% Parameters
L = 1;                % Domain length
alpha = 1;            % Thermal diffusivity (scales source term)
S0 = 1;               % Source strength
sigma = L / 40;       % Decrease Gaussian width by a factor of 10
N_max = 50;           % Maximum number of Fourier terms
x = linspace(0, L, 100); % Spatial grid

% Compute Fourier Coefficients S_n
S_n = zeros(N_max, 1);
for n = 1:N_max
    f = @(x) S0 * exp(-(x - L/2).^2 / (2 * sigma^2)) .* sin(n * pi * x / L);
    S_n(n) = (2/L) * integral(f, 0, L);
end

% Compute S(x) from Fourier series
S_reconstructed = zeros(size(x)); % Initialize reconstruction
for n = 1:N_max
    S_reconstructed = S_reconstructed + S_n(n) * sin(n * pi * x / L);
end

% Plot Original vs. Reconstructed S(x)
figure;
hold on;
plot(x, S0 * exp(-(x - L/2).^2 / (2 * sigma^2)), 'r', 'LineWidth', 2, 'DisplayName', 'Original S(x)');
plot(x, S_reconstructed, 'b--', 'LineWidth', 2, 'DisplayName', 'Reconstructed S(x)');
xlabel('x'); ylabel('S(x)');
legend;
title('Fourier Series Approximation of S(x)');
grid on;
hold off;

% Save the figure as a PNG file
saveas(gcf, 'S_reconstructed_comparison.png');

% Initialize figure for T(x) animation
fig = figure;
filename = 'heat_solution.gif';

% Compute and visualize solution convergence
T = zeros(size(x)); % Initialize temperature solution
for N = 1:N_max
    T = T + (S_n(N) / (alpha * (N * pi / L)^2)) * sin(N * pi * x / L);

    % Plotting
    plot(x, T, 'b', 'LineWidth', 2); hold on;
    plot(x, zeros(size(x)), 'k--'); % Zero reference line
    title(['Fourier Series Approximation (N = ', num2str(N), ')']);
    xlabel('x'); ylabel('Temperature T(x)');
    ylim([0 0.02]); % Adjust y-axis limits to go from 0 to 0.02
    grid on;
    hold off;

    % Capture frame for GIF
    frame = getframe(fig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if N == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end

    pause(0.05);
end

disp(['GIF saved as ', filename]);
