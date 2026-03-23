% zidane karim 3-9-2026

%% q1
x = linspace(-2*pi, 2*pi, 1000);
A = zeros(10, length(x));
current_sum = ones(size(x));
for i = 1:10
    n = 2 * i;
    current = ((-1)^(i)) * (x.^n ./ factorial(n));
    current_sum = current_sum + current;
    A(i, :) = (current_sum);
end
figure;

newColors = jet(10); % great function!
colororder(newColors);

plot(x, A);
xlabel('x');
ylabel('Approximation Value');
title('Taylor Series Approximations of cos(x) (Orders 2-20)');
grid on;
ylim([-2, 2]);
legend([arrayfun(@(i) sprintf('n = %d', 2*i), 1:10, 'UniformOutput', false), "cos(x)"]);

%% q2 
real_cos = cos(x);
figure; 

subplot(1, 2, 1);
plot(x, real_cos, 'k', 'LineWidth', 2, "Color", "#008080");
title('cos(x) Function');
grid on;

subplot(1, 2, 2);
newColors = jet(10); % great function!
colororder(newColors);
plot(x, A);
title('Taylor Series Approximations of cos(x) (Orders 2-20)');
grid on;
ylim([-2, 2]);

%% q3
S = linspace(-2*pi, 2*pi, 1000);
[X, Y] = meshgrid(S, S);

output = sinc(sqrt(X.^2+Y.^2));

figure;
surf(X, Y, output);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('sinc(sqrt(X^2 + Y^2))');
title('3D Surface Plot of sinc Function');
colorbar;
grid on;


