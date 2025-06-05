% Gauss Newton method for solving F(x) = 0 in the least squares sense
% Testing on the model f(x) = x1 * h^x2 using height/weight data

%% Define the nonlinear function and data

F = @(x, h, w) x(1) * h.^x(2) - w;

% Hardcoded Jacobian for f(x) = x1 * h^x2
DF = @(x, h) [ ...
    h.^x(2), ...
    x(1) * h.^x(2) .* log(h) ...
];

%{
% Symbolic Jacobian (optional replacement for DF)
% Requires Symbolic Math Toolbox

% syms x1 x2 h_sym
% f_sym = x1 * h_sym^x2;
% J_sym = jacobian(f_sym, [x1, x2]);
% DF = @(x, h) double( subs(J_sym, {x1, x2, h_sym}, {x(1), x(2), h(:)}) );
%}

% Height (h) and weight (w) data
h = [0.9120; 0.9860; 1.0600; 1.1300; 1.1900;
     1.2600; 1.3200; 1.3800; 1.4100; 1.4900];

w = [13.7; 15.9; 18.5; 21.3; 23.5;
     27.2; 32.7; 36.0; 38.6; 43.7];

%% Gauss Newton iteration

x = [10; 2];              % initial guess
tol = 1e-8;
max_iters = 100;
k = 0;
err = inf;

fprintf('Iter    c1              c2              ||x(k) - x(k-1)||\n');

while err > tol && k < max_iters
    k = k + 1;
    x_old = x;

    Fx = F(x, h, w);
    DFx = DF(x, h);

    lhs = DFx' * DFx;
    rhs = -DFx' * Fx;
    sk = lhs \ rhs;

    x = x + sk;
    err = norm(x - x_old);

    fprintf('%2d   %12.9f  %12.9f  %12.9f\n', k, x(1), x(2), err);
end

%% Final result

fprintf('\nFinal parameter vector\n');
fprintf('c1 = %.9f\n', x(1));
fprintf('c2 = %.9f\n', x(2));
