%% Zidane Karim -- 3/29/2026
%% q1a
n = 6;
f_lo = 10e3;
f_hi = 12e3;
r_p = 2;
f_s = 40e3;

n0 = n /2;

fcrit = [f_lo, f_hi] / (f_s / 2); % Normalize the cutoff frequencies
[z, p, k] = cheby1(n0, r_p, fcrit);
[b,a] = zp2tf(z,p,k);

figure;
title('Pole-Zero Plot', 'FontSize', 14, 'FontWeight', 'bold');
zplane(z, p); 
 
%% 1b
f = linspace(0, f_s /2, 10e4);
% normalizing the frequencies to digital radian
f_normalized = 2*pi*(f/(f_s));

H = freqz(b, a, f_normalized);          

%% 1c
% plot magnitude and phase
figure;
subplot(2,1,1);
plot(f/1e3, 20*log10(abs(H)));
xlabel('Frequency kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');
grid on;
xlim([0, (f_s/2)/1e3]);
ylim([-40, 2]);

subplot(2,1,2);
plot(f/1e3, rad2deg(unwrap(angle(H))));
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Phase Response');
xlim([0, (f_s/2)/1e3]);
grid on;



%% 1d
% using find function

magnitude_responses = 20*log10(abs(H));
lower = find(f < f_lo & magnitude_responses <= -30);
index_f1 = lower(end); 
f1 = f(index_f1);

higher = find(f > f_hi & magnitude_responses <= -30);
index_f2 = higher(1);
f2 = f(index_f2);

disp(f1);
disp(f2);

