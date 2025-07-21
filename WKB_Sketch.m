% MATLAB script to sketch triangular potential barrier
set(0,'defaultaxesfontname','Times New Roman');
set(0,'DefaultLineLineWidth', 2);
% set(0,'DefaultFigureWindowStyle','docked');
set(0,'defaultaxesfontsize',20);

clear;

% Constants
eV_to_J = 1.60218e-19; % eV to Joules

% Parameters
V0 = 10 * eV_to_J;  % Barrier height: 10 eV in Joules
E = 5 * eV_to_J;    % Particle energy: 5 eV in Joules
F = 1e-10;          % Force in N (adjusted for nm-scale width)

% Calculate turning point
x1 = (V0 - E) / F;  % Where V(x) = E

% Position range
x = linspace(-2e-9, 2 * x1, 1000); % From -2 nm to twice the barrier width
V = zeros(size(x));

% Define potential V(x)
for i = 1:length(x)
    if x(i) < 0
        V(i) = 0;       % Region I: V = 0
    else
        V(i) = V0 - F * x(i); % Region II/III: V = V0 - Fx
        if V(i) < 0
            V(i) = 0;   % Beyond barrier, V = 0
        end
    end
end

% Convert back to eV for plotting
V_eV = V / eV_to_J;
E_eV = E / eV_to_J;
V0_eV = V0 / eV_to_J;

% Plot
figure;
plot(x * 1e9, V_eV, 'b-', 'LineWidth', 2); % x in nm, V in eV
hold on;
plot([min(x) max(x)] * 1e9, [E_eV E_eV], 'r--', 'LineWidth', 1.5); % Energy line
plot([0 0], [0 V0_eV], 'k--'); % x = 0 line
plot([x1 x1] * 1e9, [0 E_eV], 'k--'); % x1 line

% Annotations
text(-1, 2, 'Region I: V = 0', 'FontSize', 14);
text(0.5 * x1 * 1e9, 8, 'Region II: Tunneling', 'FontSize', 14);
text(1.5 * x1 * 1e9, 2, 'Region III: Transmitted', 'FontSize', 14);
text(-0.5, V0_eV + 1, 'V_0 = 10 eV', 'FontSize', 14);
text(max(x) * 1e9 * 0.7, E_eV + 1, 'E = 5 eV', 'FontSize', 14);
text(x1 * 1e9, -1, 'x_1', 'FontSize', 14);
title('Triangular Potential Barrier Tunneling');
xlabel('Position x (nm)');
ylabel('Potential Energy V(x) (eV)');
grid on;
axis([min(x) * 1e9 max(x) * 1e9 -2 12]);

% Save as PNG
saveas(gcf, 'Triangular_Barrier_Sketch.png');