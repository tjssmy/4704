% MATLAB script to solve 2D harmonic wave equation in cylindrical coordinates
% Boundary condition: u(r=R, phi) = 0
% Generates geometry plot, 2x4 subplot of first 8 modes, and eigenvalue plot
% Second figure position set to [777 383 1784 855]
% Third figure plots eigenvalues

clear all;
close all;

% Parameters
R = 1; % Radius of circular domain
Nr = 30; % Number of radial intervals
Nphi = 36; % Number of azimuthal intervals
dr = R / Nr; % Radial grid spacing
dphi = 2 * pi / Nphi; % Azimuthal grid spacing
r = dr:dr:R; % Radial grid (exclude r=0 to avoid singularity)
phi = 0:dphi:(2 * pi - dphi); % Azimuthal grid for computation

% Initialize the matrix
N = (Nr - 1) * Nphi; % Number of unknowns (exclude r=R)
A = sparse(N, N);

% Build the finite difference matrix
for i = 1:Nr-1
    ri = r(i); % Current radius
    for j = 1:Nphi
        idx = (i-1) * Nphi + j; % Linear index for u(i,j)
        
        % Radial second derivative: (u_{i+1,j} - 2u_{i,j} + u_{i-1,j})/dr^2
        % Radial first derivative: (u_{i+1,j} - u_{i-1,j})/(2 dr)
        % Azimuthal second derivative: (u_{i,j+1} - 2u_{i,j} + u_{i,j-1})/(r_i^2 dphi^2)
        
        % Coefficients
        coeff_r2 = ri^2 / dr^2; % For second radial derivative
        coeff_r1 = ri / (2 * dr); % For first radial derivative
        coeff_phi = dr^2 / (ri^2 * dphi^2); % For azimuthal derivative
        
        % Diagonal term
        A(idx, idx) = -2 * coeff_r2 - 2 * coeff_phi;
        
        % Radial terms
        if i < Nr-1 % u_{i+1,j}
            idx_ip1 = i * Nphi + j;
            A(idx, idx_ip1) = coeff_r2 + coeff_r1;
        end
        if i > 1 % u_{i-1,j}
            idx_im1 = (i-2) * Nphi + j;
            A(idx, idx_im1) = coeff_r2 - coeff_r1;
        end
        
        % Azimuthal terms with periodicity
        idx_jp1 = (i-1) * Nphi + mod(j, Nphi) + 1; % u_{i,j+1}
        idx_jm1 = (i-1) * Nphi + mod(j-2, Nphi) + 1; % u_{i,j-1}
        A(idx, idx_jp1) = coeff_phi;
        A(idx, idx_jm1) = coeff_phi;
    end
end

% Boundary condition at r=R (u_{Nr,j} = 0) is implicitly handled by excluding i=Nr
% Origin (r=0) is excluded to avoid singularity

% Solve eigenvalue problem
num_modes = 8;
[V, D] = eigs(A, num_modes, 'smallestabs'); % Smallest magnitude eigenvalues
eigvals = diag(D);

% Reshape eigenvectors to 2D grid
u_modes = zeros(Nr-1, Nphi, num_modes);
for k = 1:num_modes
    u_modes(:,:,k) = reshape(V(:,k), [Nphi, Nr-1])';
end

% Create full grid for plotting (include r=R with zeros)
r_full = dr:dr:R; % Radial grid for plotting (Nr points)
u_full = zeros(Nr, Nphi, num_modes);
u_full(1:Nr-1,:,:) = u_modes; % u=0 at r=R (i=Nr)

% Extend phi for plotting to include 2*pi
phi_plot = 0:dphi:2*pi; % Nphi+1 points
u_plot = zeros(Nr, Nphi+1, num_modes);
u_plot(:,1:Nphi,:) = u_full; % Copy solution
u_plot(:,Nphi+1,:) = u_full(:,1,:); % Enforce periodicity: u(:,2*pi) = u(:,0)

% Convert to Cartesian coordinates for plotting
[Phi, R_grid] = meshgrid(phi_plot, r_full); % Create 2D grids of size Nr x (Nphi+1)
[X, Y] = pol2cart(Phi, R_grid); % Convert to Cartesian coordinates

% Plot 1: Geometry and Grid
figure(1);
% Plot boundary
theta = linspace(0, 2*pi, 100);
plot(R * cos(theta), R * sin(theta), 'b-', 'LineWidth', 2);
hold on;
% Plot grid points
plot(X(:), Y(:), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
% Plot boundary points
plot(R * cos(phi), R * sin(phi), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('x');
ylabel('y');
title('Geometry and Discretization');
legend('Boundary (u=0)', 'Grid points', 'Boundary points', 'Location', 'NorthEastOutside');
axis equal;
grid on;
hold off;
print('-dpng', 'geometry_cylindrical.png');

% Plot 2: First 8 Modes in 2x4 Subplots
figure(2);
set(gcf, 'Position', [777 383 1784 855]); % Set figure position and size
for k = 1:num_modes
    subplot(2, 4, k);
    surf(X, Y, u_plot(:,:,k), 'EdgeColor', 'none');
    colormap jet;
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Mode %d', k));
    axis equal;
    view(2);
end
sgtitle('First 8 Eigenmodes of the Harmonic Wave Equation');
print('-dpng', 'modes_cylindrical.png');

% Plot 3: Eigenvalues
figure(3);
stem(1:num_modes, abs(eigvals), 'filled', 'b', 'LineWidth', 2);
xlabel('Mode Number');
ylabel('|\lambda|');
title('Magnitudes of Eigenvalues for First 8 Modes');
grid on;
print('-dpng', 'eigenvalues_cylindrical.png');