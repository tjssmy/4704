
double_pendulum_gif

function double_pendulum_gif()
    % Parameters
    m1 = 1; m2 = 1;
    l1 = 1; l2 = 1;
    g = 9.81;

    % Initial conditions: [theta1, omega1, theta2, omega2]
    y0 = [pi/2, 0, pi/2, 0];

    % Time span
    tspan = linspace(0, 20, 1000);

    % Solve ODE
    [t, Y] = ode45(@(t, y) pendulum_ode(t, y, m1, m2, l1, l2, g), tspan, y0);

    % Prepare figure
    figure('Color', 'w');
    filename = 'double_pendulum.gif';

    % Initialize trace
    trace_x = [];
    trace_y = [];

    for i = 1:length(t)
        theta1 = Y(i,1);
        theta2 = Y(i,3);

        x1 = l1 * sin(theta1);
        y1 = -l1 * cos(theta1);
        x2 = x1 + l2 * sin(theta2);
        y2 = y1 - l2 * cos(theta2);

        % Update trace
        trace_x(end+1) = x2;
        trace_y(end+1) = y2;

        % Plot pendulum and trace
        plot([0 x1 x2], [0 y1 y2], '-o', 'LineWidth', 2, 'MarkerSize', 6);
        hold on;
        plot(trace_x, trace_y, 'r.', 'MarkerSize', 1.5);
        hold off;

        axis equal;
        axis([-2.2 2.2 -2.2 2.2]);
        title(sprintf('Double Pendulum at t = %.2f s', t(i)));
        drawnow;

        % Capture frame
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if i == 1
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
        end
    end
end

function dydt = pendulum_ode(t, y, m1, m2, l1, l2, g)
    theta1 = y(1); omega1 = y(2);
    theta2 = y(3); omega2 = y(4);

    delta = theta2 - theta1;

    den1 = (m1 + m2)*l1 - m2*l1*cos(delta)^2;
    den2 = (l2/l1)*den1;

    domega1 = (m2*l1*omega1^2*sin(delta)*cos(delta) + m2*g*sin(theta2)*cos(delta) + ...
              m2*l2*omega2^2*sin(delta) - (m1 + m2)*g*sin(theta1)) / den1;

    domega2 = (-m2*l2*omega2^2*sin(delta)*cos(delta) + (m1 + m2)*(g*sin(theta1)*cos(delta) - ...
              l1*omega1^2*sin(delta) - g*sin(theta2))) / den2;

    dydt = [omega1; domega1; omega2; domega2];
end

