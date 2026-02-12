% 1
u = -5:3:7;
v = (-pi:pi/2:pi)';
% 2
nums = 1:10;
n = prod(nums);
 % 3
A = eye(5);
A(3,3) = 0;

% B is separated by evens and odds
evens = 12:-2:2;
B_first_half = reshape(evens, 3,2); 
% not sure where i was supposed to use tranpose
B_second_half = B_first_half -1;
B = [B_first_half B_second_half];

% 4
pies = linspace(-pi, pi, 1000);
ns = (0:50)';
an = (2*ns)+1; % changed to 1-50, 0-50 results in division by zero

sins = sin(an.*pies) ./ an;
square = sum(sins);

plot(pies, square);