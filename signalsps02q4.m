% Zidane Karim 2-14-2026 ECE211
% 4c
%% figure 1
f0 = 1000; % hz
sample_rate = 5000; % hz

t_vals = 0:1/sample_rate:1;
I = 1 ./ sqrt(1+t_vals.^2);
Q = t_vals ./ sqrt(1+t_vals.^2);

x_t = I .* cos(2*pi*f0*t_vals) - Q .* sin(2*pi*f0*t_vals);

envelope = ones(size(t_vals));


figure;
plot(t_vals, x_t, 'b');
hold on;
plot(t_vals, envelope, 'r', 'LineWidth', 1.5);
hold off;

legend('x(t)', 'Envelope');
xlabel('Time (s)');
ylabel('Amplitude');
title('RF Signal and Envelope');

%% figure 2
figure;
subplot(2,1,1);
plot(t_vals, I, 'g', 'LineWidth', 1.5);
title('In-Phase Component x_I(t)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t_vals, Q, 'm', 'LineWidth', 1.5);
title('Quadrature Component x_Q(t)');
xlabel('Time (s)');
ylabel('Amplitude');