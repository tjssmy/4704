function gaussian_wavepacket()
    % Parameters
    hbar = 1; m = 1; omega = 1;
    x0 = -2; p0 = 1; sigma = 0.5;
    x = linspace(-5, 5, 200);
    t_vals = linspace(0, 10, 100);
    
    % Compute initial Gaussian wave packet
    Psi_gauss = (1/(pi * sigma^2)^(1/4)) * exp(-((x - x0).^2)/(2 * sigma^2) + 1i * p0 * x / hbar);
    
    % Compute harmonic oscillator eigenfunctions for reconstruction
    n_max = 50; % Number of states in reconstruction
    Psi_reconstruct = zeros(size(x));
    for n = 0:n_max
        Psi_n = hermiteH(n, sqrt(m * omega / hbar) * x) .* exp(-m * omega * x.^2 / (2 * hbar)) / sqrt(2^n * factorial(n));
        Psi_n = Psi_n * (m * omega / (pi * hbar))^(1/4);
        c_n = trapz(x, conj(Psi_n) .* Psi_gauss); % Compute coefficients
        Psi_reconstruct = Psi_reconstruct + c_n * Psi_n;
    end
    
    % Plot and save comparison
    figure;
    plot(x, abs(Psi_gauss).^2, 'r', 'DisplayName', 'Gaussian'); hold on;
    plot(x, abs(Psi_reconstruct).^2, 'b--', 'DisplayName', 'Reconstruction');
    xlabel('x'); ylabel('Probability Density'); legend;
    title('Gaussian Wavepacket vs. Harmonic Reconstruction');
    saveas(gcf, 'wavepacket_comparison.png');
    
    % Time evolution animation
    filename = 'wavepacket_evolution.gif';
    figure;
    for t = t_vals
        Psi_t = zeros(size(x));
        for n = 0:n_max
            Psi_n = hermiteH(n, sqrt(m * omega / hbar) * x) .* exp(-m * omega * x.^2 / (2 * hbar)) / sqrt(2^n * factorial(n));
            Psi_n = Psi_n * (m * omega / (pi * hbar))^(1/4);
            c_n = trapz(x, conj(Psi_n) .* Psi_gauss);
            Psi_t = Psi_t + c_n * Psi_n * exp(-1i * (n + 1/2) * omega * t);
        end
        
        % Plot current state
        plot(x, abs(Psi_t).^2, 'k'); hold on
        plot(x, real(Psi_t).^2, 'r');
        plot(x, imag(Psi_t).^2, 'b');

        hold off
        xlabel('x'); ylabel('Probability Density');
        title(sprintf('Time Evolution: t = %.2f', t));
        axis([-5 5 -.2 1.2]);
        pause(0.1)
        
        % Capture the frame for GIF
        frame = getframe(gcf);
        im = frame2im(frame);
        [A, map] = rgb2ind(im, 256);
            if t == t_vals(1)
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
end
