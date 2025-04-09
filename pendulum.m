% MATLAB Program to Simulate a Double Pendulum with GIF Animation and Trace

% Define physical parameters
m1 = 1;        % Mass of the first pendulum (kg)
m2 = 1;        % Mass of the second pendulum (kg)
l1 = 1;        % Length of the first pendulum (m)
l2 = 1;        % Length of the second pendulum (m)
g = 9.81;      % Gravitational acceleration (m/s^2)

% Set initial conditions
theta1 = pi/2; % Initial angle of the first pendulum (radians)
theta2 = pi/2; % Initial angle of the second pendulum (radians)
omega1 = 0;    % Initial angular velocity of the first pendulum (rad/s)
omega2 = 0;    % Initial angular velocity of the second pendulum (rad/s)

% Define time parameters
dt = 0.01;     % Time step (s)
T = 30;        % Total simulation time (s)
N = floor(T/dt); % Number of time steps

% Preallocate arrays for efficiency
t = linspace(0, T, N+1);
theta1_vec = zeros(1, N+1);
theta2_vec = zeros(1, N+1);
omega1_vec = zeros(1, N+1);
omega2_vec = zeros(1, N+1);

% Assign initial conditions
theta1_vec(1) = theta1;
theta2_vec(1) = theta2;
omega1_vec(1) = omega1;
omega2_vec(1) = omega2;

% Simulation loop using Forward Euler method
for n = 1:N
    th1 = theta1_vec(n);
    th2 = theta2_vec(n);
    om1 = omega1_vec(n);
    om2 = omega2_vec(n);
    
    alpha = th1 - th2;
    A = (m1 + m2) * l1;
    B = m2 * l2 * cos(alpha);
    D = m2 * l1 * cos(alpha);
    E = m2 * l2;
    det = A * E - B * D;
    C = m2 * l2 * om2^2 * sin(alpha) + (m1 + m2) * g * sin(th1);
    F = -m2 * l1 * om1^2 * sin(alpha) + m2 * g * sin(th2);
    theta1_dd = (-C * E + F * B) / det;
    theta2_dd = (C * D - F * A) / det;
    
    theta1_vec(n+1) = th1 + dt * om1;
    theta2_vec(n+1) = th2 + dt * om2;
    omega1_vec(n+1) = om1 + dt * theta1_dd;
    omega2_vec(n+1) = om2 + dt * theta2_dd;
end

% Compute positions for animation
x1_vec = l1 * sin(theta1_vec);l
y1_vec = -l1 * cos(theta1_vec);
x2_vec = x1_vec + l2 * sin(theta2_vec);
y2_vec = y1_vec - l2 * cos(theta2_vec);

% Animation with GIF generation and trace
gif_filename = 'double_pendulum_with_trace.gif'; % Name of the output GIF file
figure('Name', 'Double Pendulum Animation with Trace');

% Initialize arrays for trace
trace_x1 = [];
trace_y1 = [];
trace_x2 = [];
trace_y2 = [];

for n = 1:5:length(t) % Step every 5 frames to reduce GIF size
    % Append current positions to trace arrays
    trace_x1 = [trace_x1, x1_vec(n)];
    trace_y1 = [trace_y1, y1_vec(n)];
    trace_x2 = [trace_x2, x2_vec(n)];
    trace_y2 = [trace_y2, y2_vec(n)];
    
    % Plot the trace of each mass
    plot(trace_x1, trace_y1, 'r-', 'LineWidth', 1); % Trace for mass 1 (red)
    hold on;
    plot(trace_x2, trace_y2, 'b-', 'LineWidth', 1); % Trace for mass 2 (blue)
    
    % Plot the current position of the pendulum
    plot([0, x1_vec(n), x2_vec(n)], [0, y1_vec(n), y2_vec(n)], '-ko', ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    
    axis equal;
    xlim([-l1-l2-0.5, l1+l2+0.5]);
    ylim([-l1-l2-0.5, l1+l2+0.5]);
    xlabel('x (m)');
    ylabel('y (m)');
    title(['Double Pendulum Motion with Trace, t = ', num2str(t(n), '%.2f'), ' s']);
    grid on;
    drawnow;
    
    % Capture the current frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256); % Convert to indexed image with 256 colors
    
    % Write to GIF file
    if n == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.03);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
    end
    
    % Clear the current figure for the next frame
    hold off;
end