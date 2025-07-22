% pec_cylinder_scattering_TM_subplots.m
% MATLAB script to compute and plot the real parts of incident, scattered, and total electric fields (Ez)
% for a plane wave incident on a perfectly conducting cylinder (TM polarization).
% Uses eigenfunction expansion in cylindrical coordinates.
% Plots using subplots in one figure.

clear; close all; clc;

% Parameters
lambda = 1;      % Wavelength (arbitrary units)
k = 2 * pi / lambda;  % Wavenumber
b = lambda / 2;  % Cylinder radius (adjust as needed, e.g., for kb = pi)
E0 = 1;          % Incident field amplitude
Nmax = 30;       % Truncation limit for series (increase for larger kb)

% Computational grid in Cartesian coordinates
grid_size = 201; % Number of points in x and y
x = linspace(-3*b, 3*b, grid_size);
y = linspace(-3*b, 3*b, grid_size);
[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);
PHI = atan2(Y, X);

% Incident field
Ez_inc = E0 * exp(1i * k * X);

% Initialize total field
Ez_total = zeros(size(R)) + 1i * zeros(size(R));

% Compute the series sum for total field (valid for r >= b)
for n = -Nmax:Nmax
    % Precompute constants for this n
    i_n = (1i)^n;  % i^n
    Jn_kb = besselj(n, k * b);
    Hn_kb = besselh(n, 1, k * b);
    
    % Avoid division by zero or small values (though for kb > 0, Hn_kb != 0)
    if abs(Hn_kb) < eps
        continue;  % Skip if Hankel is zero (rare)
    end
    
    ratio = Jn_kb / Hn_kb;
    
    % Field contributions
    Jn_kr = besselj(n, k * R);
    Hn_kr = besselh(n, 1, k * R);
    angular_part = exp(1i * n * PHI);
    
    term = i_n * (Jn_kr - ratio * Hn_kr) .* angular_part;
    
    Ez_total = Ez_total + term;
end

Ez_total = E0 * Ez_total;

% Set field inside cylinder to zero for total
Ez_total(R < b) = 0;

% Compute scattered field
Ez_scat = Ez_total - Ez_inc;

% Set scattered field inside to zero (since not defined inside PEC)
Ez_scat(R < b) = 0;

% Create one figure with three subplots for real parts
figure('Name', 'Real Parts of Fields');
sgtitle(['Fields for kb = ' num2str(k*b)]);

subplot(1,3,1);
pcolor(X, Y, real(Ez_inc));
shading interp;
colorbar;
pbaspect([1,1,1])
% axis equal;
xlabel('x');
ylabel('y');
title('Real(E_z incident)');

subplot(1,3,2);
pcolor(X, Y, real(Ez_scat));
shading interp;
colorbar;
pbaspect([1,1,1])
% axis equal;
xlabel('x');
ylabel('y');
title('Real(E_z scattered)');

subplot(1,3,3);
pcolor(X, Y, real(Ez_total));
shading interp;
colorbar;
pbaspect([1,1,1])
% axis equal;
xlabel('x');
ylabel('y');
title('Real(E_z total)');

% Adjust figure size to make aspect ratio appropriate (wide for 1x3 subplots)
set(gcf, 'Position', [100 100 1200 300]);

% Save the figure as PNG
print('-dpng', 'fields_real_parts.png');