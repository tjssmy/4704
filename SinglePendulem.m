% Simple Pendulum Animation and GIF Generation
clear all; close all; clc;

% Parameters
g = 9.81;           % Gravitational acceleration (m/s^2)
L = 1;              % Length of the pendulum (m)
theta0 = pi/4;      % Initial angle (radians, 45 degrees)
omega0 = 0;         % Initial angular velocity (rad/s)
tspan = [0 10];     % Time span for simulation (s)
dt = 0.05;          % Time step for animation (s)

% Differential equation for the pendulum (small-angle approximation)
ode = @(t, y) [y(2); -(g/L)*sin(y(1))]; % y(1) = theta, y(2) = dtheta/dt

% Initial conditions [theta, omega]
y0 = [theta0; omega0];

% Solve the differential equation
[t, y] = ode45(ode, tspan, y0);

% Interpolate for smoother animation
t_anim = 0:dt:tspan(end);
theta = interp1(t, y(:,1), t_anim);

% Set up the figure for animation
figure;
set(gcf, 'Position', [100, 100, 600, 400]);
filename = 'pendulum_animation.gif';

% Animation loop
for i = 1:length(t_anim)
    % Calculate pendulum position
    x = L * sin(theta(i));
    y = -L * cos(theta(i));
    
    % Plot pendulum
    clf;
    hold on;
    plot([0 x], [0 y], 'b-', 'LineWidth', 2); % String
    plot(x, y, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Mass
    plot(0, 0, 'k+', 'MarkerSize', 10); % Pivot point
    axis equal;
    axis([-L-0.2 L+0.2 -L-0.2 0.2]);
    xlabel('X (m)');
    ylabel('Y (m)');
    title('Simple Pendulum Animation');
    grid on;
    
    % Capture the frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', dt);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', dt);
    end
    
    drawnow;
end

disp(['GIF saved as: ' filename]);