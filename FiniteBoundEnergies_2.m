set(0, 'defaultaxesfontname', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 20);

clear; clc;

% MATLAB program to calculate finite well bound energies

% Constants
kB = 1.380649e-23;  % Boltzmann constant (J/K)
T = 300;            % Temperature (K)
m_e = 9.1093837e-31; % Electron rest mass (kg)
m_eff = 0.067 * m_e;  % Effective mass (kg)
e = 1.60217662e-19; % Elementary charge (C)
hbar = 1.0545718e-34; % Reduced Planck's constant (J s)
sigma_v = sqrt(kB * T / m_eff);  % Standard deviation of velocity components (m/s)

Lx = 200e-9;        % Length in x (m)
Ly = 100e-9;        % Length in y (m)
lambda = 100e-9;    % Mean free path (m)

% Barrier parameters
barrier1_center = 100e-9;  % m
barrier_width = 1e-9;    % m
barrier1_left = barrier1_center - barrier_width / 2;
barrier1_right = barrier1_center + barrier_width / 2;
barrier2_left = barrier1_right + 12e-9;
barrier2_right = barrier2_left + barrier_width;
V0 = 1 * e;            % Barrier height in J
V_offset = 0.075 * e;     % Energy offset added to barrier height (J)
V0_eff = V0 + V_offset;   % Effective barrier height (J)
L_well = barrier2_left - barrier1_right;

% Function to find bound energies and plot transcendental equations
function E = find_bound_energies(m, hbar, V0, L, e, V_offset, barrier1_left, barrier1_right, barrier2_left, barrier2_right)
    a = L / 2;
    z0 = a * sqrt(2 * m * V0) / hbar;
    E = [];
    z_solutions = [];
    parity = [];
    
    % For even states
    f_even = @(z) z .* tan(z) - sqrt(z0^2 - z.^2);
    for k = 0:floor(z0 / pi - 0.5)
        z_min = k * pi;
        z_max = (k + 0.5) * pi;
        try
            z = fzero(f_even, [z_min+1e-6, z_max-1e-6]);
            en = (z / a)^2 * (hbar^2 / (2 * m));
            E = [E; en + V_offset];
            z_solutions = [z_solutions; z];
            parity = [parity; "even"];
        catch
        end
    end
    
    % For odd states
    f_odd = @(z) z .* cot(z) + sqrt(z0^2 - z.^2);
    for k = 0:floor(z0 / pi - 0.5)
        z_min = (k + 0.5) * pi;
        z_max = (k + 1) * pi;
        try
            z = fzero(f_odd, [z_min+1e-6, z_max-1e-6]);
            en = (z / a)^2 * (hbar^2 / (2 * m));
            E = [E; en + V_offset];
            z_solutions = [z_solutions; z];
            parity = [parity; "odd"];
        catch
        end
    end
    
    [E, idx] = sort(E);
    z_solutions = z_solutions(idx);
    parity = parity(idx);
    
    % Print bound energies in eV
    if ~isempty(E)
        disp('Bound energies in eV:');
        for i = 1:length(E)
            disp(E(i) / e);
        end
    else
        disp('No bound energies found.');
    end
    
    % Plotting transcendental equations
    figure;
    hold on;
    grid on;
    
    % Generate z values for plotting
    z_plot = linspace(0, z0, 1000);
    f_even_plot = zeros(size(z_plot));
    f_odd_plot = zeros(size(z_plot));
    
    for i = 1:length(z_plot)
        if z_plot(i) < z0
            f_even_plot(i) = z_plot(i) * tan(z_plot(i)) - sqrt(z0^2 - z_plot(i).^2);
            f_odd_plot(i) = z_plot(i) * cot(z_plot(i)) + sqrt(z0^2 - z_plot(i).^2);
        else
            f_even_plot(i) = NaN; % Avoid plotting beyond z0
            f_odd_plot(i) = NaN;
        end
    end
    
    % Plot even and odd functions
    plot(z_plot, f_even_plot, 'b-', 'DisplayName', 'Even: z tan(z) - sqrt(z0^2 - z^2)');
    plot(z_plot, f_odd_plot, 'r-', 'DisplayName', 'Odd: z cot(z) + sqrt(z0^2 - z^2)');
    plot([0 z0], [0 0], 'k--', 'DisplayName', 'y = 0');
    
    % Plot solutions
    for i = 1:length(z_solutions)
        z = z_solutions(i);
        if strcmp(parity(i), "even")
            y = f_even(z);
        else
            y = f_odd(z);
        end
        plot(z, y, 'ko', 'MarkerFaceColor', 'g', 'DisplayName', sprintf('Solution (%.3f eV)', E(i)/e));
        text(z, y + 0.5, sprintf('E = %.3f eV', E(i)/e), 'FontSize', 12, 'Color', 'k');
    end
    
    % Customize plot
    xlabel('z');
    ylabel('f(z)');
    title('Transcendental Equations for Finite Well Bound States');
    legend('show', 'Location', 'best');
    xlim([0 z0]);
    ylim([-10 10]); % Adjust based on expected range
    
    % Save transcendental equations plot as PNG
    saveas(gcf, 'FiniteWellPlot.png');
    
    hold off;
    
    % Plotting the potential well
    figure;
    hold on;
    grid on;
    
    % Define x values for plotting (in nm for display)
    x = linspace(barrier1_left - 10e-9, barrier2_right + 10e-9, 1000);
    V = zeros(size(x));
    
    % Define potential in joules, then convert to eV
    for i = 1:length(x)
        if x(i) < barrier1_left
            V(i) = 0;
        elseif x(i) <= barrier1_right
            V(i) = V0 + V_offset; % Corrected to include V_offset
        elseif x(i) < barrier2_left
            V(i) = V_offset;
        elseif x(i) <= barrier2_right
            V(i) = V0 + V_offset;
        else
            V(i) = 0;
        end
    end
    V_eV = V / e; % Convert to eV
    x_nm = x * 1e9; % Convert to nm for plotting
    
    % Plot potential
    plot(x_nm, V_eV, 'b-', 'DisplayName', 'Potential V(x)');
    
    % Plot bound energy levels
    for i = 1:length(E)
        energy_eV = E(i) / e;
        plot([barrier1_right * 1e9, barrier2_left * 1e9], [energy_eV, energy_eV], 'r--', 'DisplayName', sprintf('E_%d = %.3f eV', i, energy_eV));
        text(barrier2_left * 1e9 + 0.5, energy_eV, sprintf('E_%d = %.3f eV', i, energy_eV), 'FontSize', 12, 'Color', 'r');
    end
    
    % Customize plot
    xlabel('x (nm)');
    ylabel('Potential and Energy (eV)');
    title('Finite Quantum Well Potential with Bound Energy Levels');
    legend('show', 'Location', 'best');
    xlim([min(x_nm) max(x_nm)]);
    ylim([-0.05 max(V_eV) + 0.05]); % Add padding for visibility
    
    % Save potential well plot as PNG
    saveas(gcf, 'PotentialWellPlot.png');
    
    hold off;
end

% Compute bound energies
E_bound = find_bound_energies(m_eff, hbar, V0, L_well, e, V_offset, barrier1_left, barrier1_right, barrier2_left, barrier2_right);
E_bound_ev = E_bound / e;