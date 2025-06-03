% Given data
x = [0.24 0.65 0.95 1.24 1.73 2.01 2.23 2.52 2.77 2.99]';
y = [0.23 -0.26 -1.10 -0.45 0.27 0.10 -0.29 0.24 0.56 1.00]';

A = [log(x), cos(x), exp(x)];
b = y;

fprintf('\n=== Part A Output ===\n\n');

% Output the coefficient matrix and right-hand side (RHS) vector of the inconsistent linear system for a, b, and c.
fprintf('Inconsistent system: Ax = b\n\n');
for i = 1:size(A,1)
    fprintf('[%9.4f %9.4f %9.4f]', A(i,1), A(i,2), A(i,3));
    if i == size(A,1)-2
        fprintf('     [a]');
    elseif i == size(A,1)-1
        fprintf('     [b]');
    elseif i == size(A,1)
        fprintf('     [c]');
    else
        fprintf('        ');
    end
    fprintf('   =   %6.2f\n', b(i));
end

% Output the normal equations
ATA = A' * A;
ATb = A' * b;
fprintf('\nNormal equations: (A^T A)x = A^T b\n');
fprintf('[%10.4f %10.4f %10.4f]   [a]   =   [%10.4f]\n', ATA(1,1), ATA(1,2), ATA(1,3), ATb(1));
fprintf('[%10.4f %10.4f %10.4f] * [b]   =   [%10.4f]\n', ATA(2,1), ATA(2,2), ATA(2,3), ATb(2));
fprintf('[%10.4f %10.4f %10.4f]   [c]   =   [%10.4f]\n', ATA(3,1), ATA(3,2), ATA(3,3), ATb(3));

% Output the values of a, b, and c determined by solving the normal equations.
x_sol = ATA \ ATb;
fprintf('\nSolution [a; b; c] from normal equations:\n');
fprintf('a = %.10f\n', x_sol(1));
fprintf('b = %.10f\n', x_sol(2));
fprintf('c = %.10f\n', x_sol(3));

% Output the SE and RMSE for the solution from the normal equations.
y_fit = A * x_sol;
SE = sum((b - y_fit).^2);
RMSE = sqrt(SE / length(b));
fprintf('\nSquared Error (SE): %.5f\n', SE);
fprintf('Root Mean Squared Error (RMSE): %.5f\n', RMSE);

fprintf('\n=== Part B Output ===\n\n');

% (b) Determine a, b, and c using reduced QR factorization.
[m, n] = size(A);

x1 = A(:,1);
e1 = zeros(m,1); e1(1) = norm(x1)*(-sign(x1(1)));
v1 = x1 - e1;
v1 = v1 / norm(v1);
H1 = eye(m) - 2 * (v1 * v1');

A1 = H1 * A;

x2 = A1(2:end,2);
e2 = zeros(m-1,1); e2(1) = norm(x2)*(-sign(x2(1)));
v2 = x2 - e2;
v2 = v2 / norm(v2);
H2_sub = eye(m-1) - 2 * (v2 * v2');
H2 = eye(m); H2(2:end,2:end) = H2_sub;

A2 = H2 * A1;

x3 = A2(3:end,3);
e3 = zeros(m-2,1); e3(1) = norm(x3)*(-sign(x3(1)));
v3 = x3 - e3;
v3 = v3 / norm(v3);
H3_sub = eye(m-2) - 2 * (v3 * v3');
H3 = eye(m); H3(3:end,3:end) = H3_sub;

R = H3 * H2 * H1 * A;
Q = H1 * H2 * H3;

% Output the Q and R matrices obtained using Householder reflectors.
fprintf('Output the Q and R matrices obtained using Householder reflectors.\n');
fprintf('Q = \n'); disp(Q);
fprintf('R = \n'); disp(R);

x_qr = R(1:n,:) \ (Q(:,1:n)' * b);
fprintf('\nSolution [a; b; c] from reduced QR factorization:\n');
fprintf('a = %.10f\n', x_qr(1));
fprintf('b = %.10f\n', x_qr(2));
fprintf('c = %.10f\n', x_qr(3));

% Output the SE and RMSE for the solution from reduced QR factorization.
y_qr = A * x_qr;
SE_qr = sum((b - y_qr).^2);
RMSE_qr = sqrt(SE_qr / length(b));
fprintf('\nSquared Error (SE): %.5f\n', SE_qr);
fprintf('Root Mean Squared Error (RMSE): %.5f\n', RMSE_qr);

% In the same figure, plot the two fitting models along with all given data points.
xx = linspace(min(x), max(x), 300)';
A_plot = [log(xx), cos(xx), exp(xx)];
y_plot_A = A_plot * x_sol;
y_plot_B = A_plot * x_qr;

figure;
plot(x, b, 'ko', 'MarkerFaceColor', 'k'); hold on;
plot(xx, y_plot_A, 'b-', 'LineWidth', 1.5);
plot(xx, y_plot_B, 'r--', 'LineWidth', 1.5);
legend('Data Points', 'Normal Equation Fit', 'QR Fit', 'Location', 'Best');
xlabel('x'); ylabel('y');
title('Least Squares Fits: Normal Equations vs QR Factorization');
grid on;
