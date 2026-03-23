% Zidane Karim -- 3/23/2026
%% q1

fs = 44100; % sample
n = 88200;
T = 1 / fs;
t = (0:n-1) * T;

freqs = [14.96e3, 5.52e3, 3.15e3, -3.76e3, 11.01e3];
dbm_freqs = [20, 19, 7, -7, -4];

% convertinjg to deicbels
scale = 10.^(dbm_freqs / 20);


x = zeros(1, n);

for i = 1:length(freqs)
    x = x + dbm_freqs(i) * exp(1j * 2 * pi * freqs(i) * t);
end

noise_dBm = -10;
noise_mag = 10^(noise_dBm / 20);
noise = noise_mag * randn(1, n); % random magnitude now
x = x + noise; 

X = fft(x);
X_shifted = fftshift(X);

frequency_axis = linspace(-fs/2, fs/2 - fs/n, n);
X_mag = abs(X_shifted) / n;
X_dB = 20 * log10(X_mag);

figure('Name', 'Audio Sample');

subplot(3, 1, 1);
plot(t, real(x));
title('Real Component of Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.005]); 

subplot(3, 1, 2);
plot(t, imag(x));
title('Imaginary Component of Signal');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0, 0.005]);

subplot(3, 1, 3);
plot(frequency_axis, X_dB);
title('DFT of Signal (Decibel Scale)');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

