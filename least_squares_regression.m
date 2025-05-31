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
