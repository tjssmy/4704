% MATLAB program to create, analyze, and reconstruct a modulated Gaussian electron wave function
% with spatial offset and animation

clear all;
close all;

% Parameters
L = 100;              % Spatial domain length (arbitrary units)
N = 1024;            % Number of points (power of 2 for FFT)
dx = L/N;            % Spatial step size
x = 0:dx:L-dx;  % Spatial coordinate array

% Wave function parameters
sigma = 5;           % Gaussian width (spatial spread)
k0 = 0.5;           % Central wave number (modulation frequency)
A = 1;              % Amplitude
x0 = 30;           % Spatial offset

% Create the modulated Gaussian pulse (wave function) with offset
psi = A * exp(-((x-x0).^2)/(2*sigma^2)) .* exp(1i*k0*(x-x0));
% psi = cos(2*pi*x/L*20);

% Calculate the FFT for spectral content
psi_k = fftshift(fft(psi));
dk = 2*pi/L;         % Frequency step size
k = -N/2*dk:dk:(N/2-1)*dk;  % Wave number array



% Calculate magnitudes for plotting
psi_mag = abs(psi);         % Magnitude of spatial wave function
psi_k_mag = abs(psi_k);     % Magnitude of spectral content

% Normalize spectral content for better visualization
psi_k_mag = psi_k_mag/max(psi_k_mag);

% Create figure with two subplots
figure('Position', [100 100 800 400]);

% Plot spatial domain (wave function)
subplot(1,2,1);
plot(x, psi_mag, 'b-', 'LineWidth', 1.5);
hold on;
plot(x, real(psi), 'r--', 'LineWidth', 1);
plot(x, imag(psi), 'g--', 'LineWidth', 1);
hold off;
grid on;
xlabel('Position (x)');
ylabel('Amplitude');
title('Electron Wave Function (Offset by -30)');
legend('Magnitude', 'Real', 'Imaginary');
axis tight;

% Plot frequency domain (spectral content)
subplot(1,2,2);
plot(k, psi_k_mag, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Wave number (k)');
ylabel('Normalized Amplitude');
title('Spectral Content');
axis tight;

% Adjust figure properties
set(gcf, 'Color', 'white');
sgtitle('Modulated Gaussian Electron Wave Function Analysis (Offset)');

% Display some basic properties
disp('Wave function properties:');
disp(['Spatial width (sigma): ', num2str(sigma)]);
disp(['Central wave number (k0): ', num2str(k0)]);
disp(['Spatial offset (x0): ', num2str(x0)]);
disp(['Uncertainty product (Δx·Δk): ', num2str(sigma * (1/sigma))]);

% Reconstruction using inverse Fourier synthesis
psi_reconstructed = zeros(size(x));  % Initialize reconstructed wave function
n_steps = 50;                       % Number of steps for animation
step_size = floor(N/n_steps);       % Number of components per step

% Prepare animation figure
figure('Position', [100 100 600 400]);
h_real = plot(x, zeros(size(x)), 'r-', 'LineWidth', 1.5);
hold on;
h_imag = plot(x, zeros(size(x)), 'g-', 'LineWidth', 1.5);
plot(x, real(psi), 'r--', 'LineWidth', 1);
plot(x, imag(psi), 'g--', 'LineWidth', 1);
hold off;
grid on;
xlabel('Position (x)');
ylabel('Amplitude');
title('Wave Function Reconstruction');
legend('Reconstructed Real', 'Reconstructed Imag', 'Original Real', 'Original Imag');
axis([min(x) max(x) -1.5*A 1.5*A]);
set(gcf, 'Color', 'white');

% Animation and GIF creation
filename = 'wave_reconstruction.gif';
for n = 1:n_steps
    % Add next batch of spectral components
    start_idx = 1 + (n-1)*step_size;
    end_idx = min(n*step_size, N);
    
    for m = start_idx:end_idx
        % Correct k-index mapping
        k_idx = m - 1 - N/2;  % Proper centering of k values
        psi_reconstructed = psi_reconstructed + ...
            (psi_k(m)/N) * exp(1i * k(k_idx + N/2 + 1) * x);
    end
    
    % Update plot with real and imaginary parts
    set(h_real, 'YData', real(psi_reconstructed));
    set(h_imag, 'YData', imag(psi_reconstructed));
    drawnow;
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF file
    if n == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

disp('Animation complete. GIF file saved as wave_reconstruction.gif');

% Additional parameters for proper phase propagation
t_max = 50;          % Total time for animation
n_frames = 100;      % Number of frames
dt = t_max/n_frames; % Time step
m = 1;              % Mass (arbitrary units)
hbar = 1;           % Reduced Planck's constant (arbitrary units)

% Store original spectral components
psi_k_original = fftshift(fft(psi));  % From original code

% Prepare figure for propagation animation
figure('Position', [100 100 600 400]);
h_mag = plot(x, zeros(size(x)), 'b-', 'LineWidth', 1.5);
hold on;
h_real = plot(x, zeros(size(x)), 'r--', 'LineWidth', 1);
h_imag = plot(x, zeros(size(x)), 'g--', 'LineWidth', 1);
hold off;
grid on;
xlabel('Position (x)');
ylabel('Amplitude');
title('Propagating Gaussian Wave Packet');
legend('Magnitude', 'Real', 'Imaginary');
axis([min(x) max(x) -1.5*A 1.5*A]);
set(gcf, 'Color', 'white');

% Animation of propagating pulse with proper phase velocities
filename_prop = 'wave_propagation.gif';

for n = 1:n_frames
    % Current time
    t = (n-1)*dt;
    
    % Initialize time-evolved spectral components
    psi_k_t = psi_k_original;
    
    % Apply phase evolution to each k-component
    for j = 1:N
        % Get corresponding wave number
        k_idx = j - 1 - N/2;  % Center k array
        k_val = k(k_idx + N/2 + 1);
        
        % Dispersion relation: ω(k) = ħk²/2m
        omega = (hbar * k_val^2) / (2*m);
        
        % Phase evolution: exp(-iωt)
        psi_k_t(j) = psi_k_original(j) * exp(-1i * omega * t);
    end
    
    % Inverse FFT to get time-evolved spatial wave function
    psi_t = ifft(ifftshift(psi_k_t));
    
    % Update plot data
    set(h_mag, 'YData', abs(psi_t));
    set(h_real, 'YData', real(psi_t));
    set(h_imag, 'YData', imag(psi_t));
    
    % Add time stamp to title
    title(['Propagating Gaussian Wave Packet (t = ' num2str(t, '%.1f') ')']);
    drawnow;
    
    % Capture frame for GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to GIF file
    if n == 1
        imwrite(imind, cm, filename_prop, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, filename_prop, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

disp('Propagation animation complete. GIF file saved as wave_propagation.gif');

% Calculate and display group velocity
k_center = k0;  % Central wave number
vg = (hbar * k_center) / m;  % Group velocity: dω/dk = ħk/m
disp(['Calculated group velocity: ' num2str(vg)]);