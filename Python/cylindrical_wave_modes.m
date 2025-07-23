function cylindrical_wave_modes()
% CYLINDRICAL_WAVE_MODES - Solve 2D harmonic wave equation in cylindrical coordinates
%
% Solves: ∇²u = -k²u in cylindrical domain with u(R,φ) = 0
% Using separation of variables and Bessel functions
%
% Generates:
% 1. Geometry and discretization figure
% 2. First 8 eigenmodes in 2x4 subplot arrangement

%% Parameters
R = 1.0;           % Domain radius
nr = 50;           % Radial grid points
nphi = 100;        % Angular grid points for plotting

%% Create coordinate grids
r = linspace(0, R, nr);
phi = linspace(0, 2*pi, nphi);
[R_mesh, PHI_mesh] = meshgrid(r, phi);
X_mesh = R_mesh .* cos(PHI_mesh);
Y_mesh = R_mesh .* sin(PHI_mesh);

%% Figure 1: Geometry and Grid
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Domain geometry
subplot(1, 2, 1);
theta = linspace(0, 2*pi, 200);
x_boundary = cos(theta);
y_boundary = sin(theta);

plot(x_boundary, y_boundary, 'k-', 'LineWidth', 3);
hold on;
fill(x_boundary, y_boundary, [0.9, 0.95, 1], 'FaceAlpha', 0.3);

% Add coordinate system
arrow_props = {'LineWidth', 2, 'Color', 'red', 'MaxHeadSize', 0.8};
quiver(0, 0, 0.7, 0, 0, arrow_props{:});
quiver(0, 0, 0, 0.7, 0, arrow_props{:});
text(0.75, -0.08, 'x', 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');
text(-0.08, 0.75, 'y', 'FontSize', 14, 'Color', 'red', 'FontWeight', 'bold');

% Add labels and annotations
text(0, -0.3, 'Ω', 'FontSize', 16, 'HorizontalAlignment', 'center');
text(0.7, 0.7, 'u(R,φ) = 0', 'FontSize', 12, 'BackgroundColor', 'white', ...
     'EdgeColor', 'black', 'Margin', 3);

xlim([-1.3, 1.3]);
ylim([-1.3, 1.3]);
axis equal;
grid on;
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
title({'Cylindrical Domain', '∇²u = -k²u'}, 'FontSize', 14, 'FontWeight', 'bold');

% Subplot 2: Grid discretization
subplot(1, 2, 2);

% Plot radial grid lines
nr_plot = 8;  % Fewer lines for clarity
nphi_plot = 16;
r_plot = linspace(0, R, nr_plot);
phi_plot = linspace(0, 2*pi, nphi_plot+1);

for i = 2:nr_plot  % Skip r=0
    theta_circle = linspace(0, 2*pi, 100);
    x_circle = r_plot(i) * cos(theta_circle);
    y_circle = r_plot(i) * sin(theta_circle);
    plot(x_circle, y_circle, 'g-', 'LineWidth', 0.8, 'Color', [0, 0.7, 0]);
    hold on;
end

% Plot angular grid lines
for j = 1:nphi_plot
    r_line = linspace(0, R, 20);
    x_line = r_line * cos(phi_plot(j));
    y_line = r_line * sin(phi_plot(j));
    plot(x_line, y_line, 'b-', 'LineWidth', 0.8, 'Color', [0, 0, 0.7]);
end

% Boundary
plot(x_boundary, y_boundary, 'k-', 'LineWidth', 3);

% Grid points
[R_grid, PHI_grid] = meshgrid(r_plot(2:end), phi_plot(1:end-1));
X_grid = R_grid .* cos(PHI_grid);
Y_grid = R_grid .* sin(PHI_grid);
plot(X_grid(:), Y_grid(:), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'red');

xlim([-1.3, 1.3]);
ylim([-1.3, 1.3]);
axis equal;
grid on;
xlabel('x', 'FontSize', 12);
ylabel('y', 'FontSize', 12);
title({'Finite Difference Grid', '(r,φ) coordinates'}, 'FontSize', 14, 'FontWeight', 'bold');

% Add legend
legend({'Circular grid (r)', 'Radial grid (φ)', 'Boundary', 'Grid points'}, ...
       'Location', 'northwest', 'FontSize', 10);

sgtitle('2D Cylindrical Wave Equation: Geometry and Discretization', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'geometry_discretization.png');

%% Calculate analytical eigenvalues and modes
% First few zeros of Bessel functions
J0_zeros = [2.4048, 5.5201, 8.6537, 11.7915];  % J₀ zeros
J1_zeros = [3.8317, 7.0156, 10.1735, 13.3237]; % J₁ zeros  
J2_zeros = [5.1356, 8.4172, 11.6198, 14.7960]; % J₂ zeros
J3_zeros = [6.3802, 9.7610, 13.0152, 16.2235]; % J₃ zeros

% Collect mode data
modes = [];

% (0,1), (0,2), (0,3) - Axisymmetric modes
for n = 1:3
    k_val = J0_zeros(n) / R;
    mode_data.m = 0;
    mode_data.n = n;
    mode_data.k = k_val;
    mode_data.type = 'axisymmetric';
    modes = [modes, mode_data];
end

% (1,1), (1,2), (1,3) - First azimuthal modes  
for n = 1:3
    k_val = J1_zeros(n) / R;
    mode_data.m = 1;
    mode_data.n = n;
    mode_data.k = k_val;
    mode_data.type = 'cos';
    modes = [modes, mode_data];
end

% (2,1), (2,2) - Second azimuthal modes
for n = 1:2
    k_val = J2_zeros(n) / R;
    mode_data.m = 2;
    mode_data.n = n;
    mode_data.k = k_val;
    mode_data.type = 'cos';
    modes = [modes, mode_data];
end

%% Figure 2: Mode shapes (2x4 subplot)
figure('Position', [150, 150, 1600, 800]);

for idx = 1:8
    subplot(2, 4, idx);
    
    mode = modes(idx);
    m = mode.m;
    n = mode.n;
    k = mode.k;
    
    % Calculate mode shape
    U = zeros(size(R_mesh));
    
    for i = 1:length(r)
        for j = 1:length(phi)
            r_val = r(i);
            phi_val = phi(j);
            
            if r_val == 0
                if m == 0
                    radial_part = 1;  % J₀(0) = 1
                else
                    radial_part = 0;  % Jₘ(0) = 0 for m > 0
                end
            else
                radial_part = besselj(m, k * r_val);
            end
            
            if strcmp(mode.type, 'axisymmetric')
                angular_part = 1;
            else
                angular_part = cos(m * phi_val);
            end
            
            U(j, i) = radial_part * angular_part;
        end
    end
    
    % Plot
    contourf(X_mesh, Y_mesh, U, 20, 'LineStyle', 'none');
    hold on;
    contour(X_mesh, Y_mesh, U, 10, 'k-', 'LineWidth', 0.5, 'LineStyle', '-');
    
    % Add boundary circle
    plot(x_boundary, y_boundary, 'k-', 'LineWidth', 2);
    
    % Formatting
    axis equal;
    xlim([-1.1, 1.1]);
    ylim([-1.1, 1.1]);
    
    % Colormap and colorbar
    colormap(redblue_colormap());
    c = colorbar;
    c.Label.String = 'u(r,φ)';
    c.Label.FontSize = 10;
    
    % Title
    if m == 0
        title_str = sprintf('Mode (0,%d)\nk = %.3f', n, k);
    else
        title_str = sprintf('Mode (%d,%d)\nk = %.3f', m, n, k);
    end
    title(title_str, 'FontSize', 11, 'FontWeight', 'bold');
    
    % Remove ticks for cleaner look
    set(gca, 'XTick', [], 'YTick', []);
    
    % Add mode pattern description
    if m == 0
        pattern_desc = 'Axisymmetric';
    elseif m == 1
        pattern_desc = 'Dipole';
    elseif m == 2
        pattern_desc = 'Quadrupole';
    else
        pattern_desc = sprintf('%d-pole', 2*m);
    end
    
    text(0, -1.25, pattern_desc, 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.5, 0.5, 0.5]);
end

sgtitle('Cylindrical Wave Equation Eigenmodes: u(r,φ) = J_m(k_{mn}r)[cos(mφ)/sin(mφ)]', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'modes_cylindrical.png');

%% Figure 3: Eigenvalue Analysis
figure('Position', [200, 200, 1400, 600]);

% Extend mode data for better visualization
all_modes = [];

% Calculate more modes for each m
m_values = [0, 1, 2, 3, 4];
n_max = 6;  % Number of radial modes to calculate

for m_idx = 1:length(m_values)
    m = m_values(m_idx);
    
    % Get Bessel function zeros
    if m == 0
        zeros_list = [2.4048, 5.5201, 8.6537, 11.7915, 15.2049, 18.4445];
    elseif m == 1
        zeros_list = [3.8317, 7.0156, 10.1735, 13.3237, 16.4599, 19.5808];
    elseif m == 2
        zeros_list = [5.1356, 8.4172, 11.6198, 14.7960, 17.9597, 21.1169];
    elseif m == 3
        zeros_list = [6.3802, 9.7610, 13.0152, 16.2235, 19.4148, 22.5982];
    elseif m == 4
        zeros_list = [7.5883, 11.0647, 14.3725, 17.6160, 20.8269, 24.0182];
    end
    
    for n = 1:n_max
        mode_info.m = m;
        mode_info.n = n;
        mode_info.k = zeros_list(n) / R;
        all_modes = [all_modes, mode_info];
    end
end

% Subplot 1: Eigenvalues vs radial mode number for different m
subplot(1, 3, 1);
colors = {'b-o', 'r-s', 'g-^', 'm-d', 'c-v'};
hold on;

for m_idx = 1:length(m_values)
    m = m_values(m_idx);
    n_vals = [];
    k_vals = [];
    
    for mode = all_modes
        if mode.m == m
            n_vals = [n_vals, mode.n];
            k_vals = [k_vals, mode.k];
        end
    end
    
    plot(n_vals, k_vals, colors{m_idx}, 'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', sprintf('m = %d', m));
end

xlabel('Radial mode number (n)', 'FontSize', 12);
ylabel('Eigenvalue k_{mn}', 'FontSize', 12);
title('Eigenvalues vs Radial Mode Number', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on;
xlim([0.5, n_max + 0.5]);

% Subplot 2: Eigenvalues vs azimuthal mode number for different n
subplot(1, 3, 2);
colors_n = {'b-o', 'r-s', 'g-^', 'm-d', 'c-v', 'k-p'};
hold on;

for n = 1:4  % First 4 radial modes
    m_vals = [];
    k_vals = [];
    
    for mode = all_modes
        if mode.n == n && mode.m <= 4
            m_vals = [m_vals, mode.m];
            k_vals = [k_vals, mode.k];
        end
    end
    
    plot(m_vals, k_vals, colors_n{n}, 'LineWidth', 2, 'MarkerSize', 6, ...
         'DisplayName', sprintf('n = %d', n));
end

xlabel('Azimuthal mode number (m)', 'FontSize', 12);
ylabel('Eigenvalue k_{mn}', 'FontSize', 12);
title('Eigenvalues vs Azimuthal Mode Number', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
grid on;
xlim([-0.5, 4.5]);

% Subplot 3: 2D eigenvalue map
subplot(1, 3, 3);

% Create matrices for surface plot
m_range = 0:4;
n_range = 1:6;
[M_grid, N_grid] = meshgrid(m_range, n_range);
K_grid = zeros(size(M_grid));

for i = 1:length(n_range)
    for j = 1:length(m_range)
        m = M_grid(i, j);
        n = N_grid(i, j);
        
        % Find corresponding eigenvalue
        for mode = all_modes
            if mode.m == m && mode.n == n
                K_grid(i, j) = mode.k;
                break;
            end
        end
    end
end

% Create surface plot
surf(M_grid, N_grid, K_grid, 'EdgeColor', 'black', 'LineWidth', 0.5);
colormap(jet);
colorbar;

% Add data points
hold on;
for mode = all_modes
    if mode.m <= 4 && mode.n <= 6
        plot3(mode.m, mode.n, mode.k, 'ko', 'MarkerSize', 8, ...
              'MarkerFaceColor', 'white', 'LineWidth', 2);
    end
end

xlabel('Azimuthal mode (m)', 'FontSize', 12);
ylabel('Radial mode (n)', 'FontSize', 12);
zlabel('Eigenvalue k_{mn}', 'FontSize', 12);
title('Eigenvalue Surface k_{mn}', 'FontSize', 14, 'FontWeight', 'bold');
view(45, 30);
grid on;

% Add text annotations for some key eigenvalues
text(0, 1, all_modes(1).k + 1, sprintf('k_{01} = %.3f', all_modes(1).k), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');
text(1, 1, all_modes(7).k + 1, sprintf('k_{11} = %.3f', all_modes(7).k), ...
     'FontSize', 10, 'FontWeight', 'bold', 'Color', 'red');

sgtitle('Cylindrical Wave Equation: Eigenvalue Analysis', ...
        'FontSize', 16, 'FontWeight', 'bold');

% Save figure
saveas(gcf, 'eigenvalues_cylindrical.png');

%% Display eigenvalue table
fprintf('\n=== Cylindrical Wave Equation Eigenvalues ===\n');
fprintf('Mode (m,n)    k_analytical    Pattern\n');
fprintf('----------------------------------------\n');
for idx = 1:8
    mode = modes(idx);
    if mode.m == 0
        fprintf('  (0,%d)         %.6f      Axisymmetric\n', mode.n, mode.k);
    else
        pattern_name = {'Dipole', 'Quadrupole'};
        if mode.m <= 2
            pattern = pattern_name{mode.m};
        else
            pattern = sprintf('%d-pole', 2*mode.m);
        end
        fprintf('  (%d,%d)         %.6f      %s\n', mode.m, mode.n, mode.k, pattern);
    end
end
fprintf('----------------------------------------\n');

fprintf('\nGenerated files:\n');
fprintf('- geometry_discretization.png\n');
fprintf('- modes_cylindrical.png\n');
fprintf('- eigenvalues_cylindrical.png\n');

end

function cmap = redblue_colormap()
% Create a red-white-blue colormap
n = 64;
r = [linspace(0, 1, n/2), ones(1, n/2)]';
g = [linspace(0, 1, n/2), linspace(1, 0, n/2)]';
b = [ones(1, n/2), linspace(1, 0, n/2)]';
cmap = [r, g, b];
end
