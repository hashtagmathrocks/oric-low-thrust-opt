clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indirect method
% LEO: omega = 4 rad/h
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global omega alpha rho A B x0 xf t0 tf
% Constant
omega = 4;                                  % angular velocity, 4 rad/h
alpha = 1e5;                                % Parameter to be adjusted
rho = 10;                                   % Distance between chief and deputy

% Initial and final states
x0 = [-rho; 0; 0; 0; 0; pi];
xf = [0; -rho; 0; 0; 0; pi];
t0 = 0;
tf = 0.25;

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
A = [zeros(3), eye(3); M1, M2];
B = [zeros(3); eye(3)];

lambda_guess = [7.79e4; 9.37e4; 3.35e4; -456.09; 384.13; 165.61];

% x_sol is lambda
x_sol = fsolve(@solver, lambda_guess);

lambda = x_sol;

[t, y] = ode45(@bvpfun, [t0, tf], [x0', lambda']);
lambda = y(:, 7:12)';
u = B' * lambda;
dJ = 0.5 .* (vecnorm(u) .* vecnorm(u));
J = trapz(t, dJ);
fprintf('J = %f', J);

% Plot
subplot(2,2,1)
plot(t, x, 'LineWidth', 1.5);
title('State - pos');

% Costate
costate = lambda;
subplot(2,2,2)
plot(t, costate(1, :), 'LineWidth', 1.5);hold on
plot(t, costate(2, :), 'LineWidth', 1.5);hold on
plot(t, costate(3, :), 'LineWidth', 1.5);hold on
plot(t, costate(4, :), 'LineWidth', 1.5);hold on
plot(t, costate(5, :), 'LineWidth', 1.5);hold on
plot(t, costate(6, :), 'LineWidth', 1.5);
legend('costate1', 'costate2', 'costate3', 'costate4', 'costate5', 'costate6');
title('Costate');

% Control
subplot(2,2,3)
plot(t, dJ, 'LineWidth', 1.5);
title('Control');

% Trajectory
subplot(2,2,4)
r = 10;
[X, Y, Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'bone')
axis equal
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Chief');hold on
plot3(x1(1), x2(1), x3(1), 'g*', 'LineWidth', 2);hold on
text(x1(1), x2(1), x3(1), 'Departure');hold on
plot3(x1(end), x2(end), x3(end), 'c*', 'LineWidth', 2);hold on
text(x1(end), x2(end), x3(end), 'Arrival');hold on
plot3(x1, x2, x3, 'k-', 'LineWidth', 1.5);hold on
quiver3(x1, x2, x3, u(1, :), u(2, :), u(3, :), 0.3, 'Color', 'r','LineWidth', 1);
title('Trajectory');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control equations
% f = dydt
% y = [x; lambda], 12x1 vector
% x: state - v and a, 6x1 vector
% lambda: costate, 6x1 vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = solver(x)
global x0 xf t0 tf
initial_guess = [x0', x'];
options = odeset('MaxStep', 1e-4, 'RelTol',1e-5,'AbsTol',1e-9);
[t, y] = ode45(@bvpfun, [t0, tf], initial_guess);
s = length(t);
F = y(s, 1:6) - xf';
end

function dydt = bvpfun(t, y)
global rho alpha A B
dydt = zeros(12, 1);
syms x1 x2 x3 p1 p2 p3
C =x1^2 + x2^2 + x3^2 - rho^2;
C_value = double(subs(C, {x1, x2, x3}, {y(1), y(2), y(3)}));
gradC = [diff(C, x1); diff(C, x2); diff(C, x3); diff(C, p1); diff(C, p2); diff(C, p3)];
gradC_value = double(subs(gradC, {x1, x2, x3, p1, p2, p3}, {y(1), y(2), y(3), y(4), y(5), y(6)}));
dydt(1:6) = A * y(1:6) - B * B' * y(7:12);
dydt(7:12) = -A' * y(7:12) - alpha * C_value * gradC_value;
end