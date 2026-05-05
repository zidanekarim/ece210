% Zidane Karim PS08 4/30

%% Part 1
N = 1e6;
x_gauss = randn(N, 1);

x_student_t = trnd(5, N, 1) * sqrt(3/5); % sqrt(3/5) comes from for a t distribution with nu>2, variance = nu/(nu-2). Here nu = 5 -> var = 5/(5-2) = 5/3.
% To convert to variance 1 multiply by sqrt(1 / (5/3)) = sqrt(3/5).

% cauchy factor given
U = rand(N, 1);
x_cauchy = 0.544 * tan(pi*U);

frac_gauss = sum(abs(x_gauss) < 1) / N;
frac_t = sum(abs(x_student_t) < 1) / N;
frac_cauchy = sum(abs(x_cauchy) < 1) / N;

fprintf('Fraction |X| < 1:\nGaussian: %.4f\nStudent''s t: %.4f\nCauchy: %.4f\n', ...
    frac_gauss, frac_t, frac_cauchy);

figure;
subplot(3,1,1); plot(x_gauss); title('Gaussian N(0,1)');
yline(1, '--r'); yline(-1, '--r');
subplot(3,1,2); plot(x_student_t); title('Student''s t (v=5, scaled)');
yline(1, '--r'); yline(-1, '--r');
subplot(3,1,3); plot(x_cauchy); title('Cauchy (\alpha=0.544)');
yline(1, '--r'); yline(-1, '--r');

%  b
segments_gauss = reshape(x_gauss, 100000, 10);
segments_t = reshape(x_student_t, 100000, 10);
segments_cauchy = reshape(x_cauchy, 100000, 10);

mean_gauss = mean(segments_gauss);
mean_t = mean(segments_t);
mean_cauchy = mean(segments_cauchy);

disp('Segment Means for Cauchy:');
disp(mean_cauchy);

% NOT Valid to say it is converging to 0 since there is random oscillation

%% Part 2
% a) poles and zeroes
MA = [1, 0.4, 0.2];
AR = [1, -1.6, 0.81];

zeros_ROOTS = roots(MA);
poles_ROOTS = roots(AR);
figure; 
zplane(MA, AR); 
title('Pole-Zero Plot of H(z)');
% b
N_samples=1e4;
sigma_v = sqrt(2);
v = sigma_v * randn(N_samples, 1); 
x = filter(MA, AR, v); % Apply innovations filter

r_x = zeros(7, 1);
for m = 0:6
    r_x(m+1) = dot(x(1:N_samples-m), x(m+1:N_samples)) / N_samples;
end

r_x_full = [flipud(r_x(2:end)); r_x]; % flipud reflects about x-axis, needed bc symmetry
m_vec = -6:6;

figure; 
stem(m_vec, r_x_full, 'filled'); 
title('Estimated Autocorrelation r_x[m]'); 
xlabel('Lag of m');

R1 = toeplitz(r_x);
disp('Toeplitz Matrix R1:'); 
disp(R1);

eigenvalues_R1 = eig(R1);
if all(eigenvalues_R1 > 0)
    disp('R1 is Positive Definite (all eigenvalues > 0).');
end

X_mat = zeros(N_samples - 6, 7); % N - 6 because we want consistency btwn lags
for k = 1:7
    X_mat(:, k) = x(7-k+1 : N_samples-k+1);
end
R2 = (X_mat' * X_mat) / (N_samples - 6);
% this is b.6
matrix_diff_norm = norm(R1 - R2);
fprintf('Difference norm ||R1 - R2||: %f\n', matrix_diff_norm);

% part c
[s_est, w] = pwelch(x, hamming(512), 256, 512);
figure; 
plot(w, s_est);
title('Estimated PSD S_x(\omega) using Welch '); 
xlabel('Linear Normalized Frequency (\times \pi rad/sample)'); ylabel('Power');

[placeholder, peak_idx] = max(s_est);
w0 = w(peak_idx);
pole_angle = abs(angle(poles_ROOTS(1)));

fprintf('Rough Peak Frequency w0: %f rad/sample\n', w0);
fprintf('Pole Angle: %f rad/sample\n', pole_angle);

% part d
[a, varv] = aryule(x,4);
fprintf('Compare varv: %f (to true variance = 2)\n', varv);

x_ar = filter(1, a, v);

% -6 <= m <= 6
r_ar = zeros(7, 1);
for m = 0:6
    r_ar(m+1) = dot(x_ar(1:N_samples-m), x_ar(m+1:N_samples)) / N_samples;
end
r_ar_full = [flipud(r_ar(2:end)); r_ar]; % again we flip

figure;
stem(m_vec, r_x_full, 'filled', 'Color', 'b');
hold on;
stem(m_vec, r_ar_full, 'Marker', 'x', 'Color', 'r');
title('Correlation Comparison: True ARMA vs AR(4) Model');
legend('Original Signal', 'AR(4) Model'); xlabel('Lag m');
disp("At m>3, we can start to see the two graphs grow apart");
figure;
% superimposed stem plot for x0 and x
stem(1:100, x(1:100), 'filled', 'Color', 'b');
hold on;
stem(1:100, x_ar(1:100), 'Marker', 'x', 'Color', 'r');
title('First 100 Data Points Comparison');
legend('Original x[n]', 'AR(4) generated x[n]');
xlabel('Time Index n');

disp("Yes, first 100 points dont match but seem similar");