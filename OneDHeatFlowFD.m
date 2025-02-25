% Parameters
L = 1;                   % Length of the domain
alpha = 1;               % Thermal diffusivity
S0 = 1;                  % Source strength
sigma = L / 40;          % Width of the Gaussian source
N_x = 100;               % Number of spatial grid points
dx = L / (N_x - 1);      % Spatial step size
N_t = 6000;              % Increased time steps (3x longer simulation)
S_max = 1;               % Maximum source strength
x = linspace(0, L, N_x); % Spatial grid

% Adjust dt to ensure stability (CFL condition)
dt = (dx^2) / (2 * alpha);  % Stability condition for explicit method

% Compute source term S(x) as Gaussian
S = S_max * exp(-(x - L/2).^2 / (2 * sigma^2));

% Precompute the constant for the finite difference scheme
r = alpha * dt / dx^2;

% Initialize the temperature field (initial condition)
T = zeros(N_x, 1);

% Create a figure for visualization
fig = figure;

% Initialize the filename for the GIF
gif_filename = 'heat_solution_fd.gif';

% Time-stepping loop for the finite difference solution
for n = 1:N_t
    T_new = T;  % Copy of the current temperature for updates

    % Update the interior points using finite difference equation
    for i = 2:N_x-1
        T_new(i) = T(i) + r * (T(i+1) - 2*T(i) + T(i-1)) + dt * S(i);
    end

    % Apply boundary conditions (T(0,t) = T(L,t) = 0)
    T_new(1) = 0;  % Left boundary
    T_new(N_x) = 0; % Right boundary

    % Update the temperature field
    T = T_new;

    % Plotting and saving frames for GIF (only every 30th time step)
    if mod(n, 30) == 0  % Save every 30th frame
        plot(x, T, 'b', 'LineWidth', 2);
        title(['Temperature Distribution at Time Step ', num2str(n)]);
        xlabel('x');
        ylabel('Temperature T(x,t)');
        ylim([0 0.02]); % Set y-axis limits for better visualization
        grid on;

        % Capture frame for GIF
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        % Write the first frame to initialize the GIF
        if n == 30  % Start appending after the first frame
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            % Append subsequent frames to the GIF
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end

        pause(0.05); % Pause to visualize the plot
    end
end

disp(['GIF saved as ', gif_filename]);
