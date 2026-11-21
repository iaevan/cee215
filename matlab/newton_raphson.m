function [root, iter, history, flag] = newton_raphson(f, df, x0, es, maxit)
% NEWTON_RAPHSON  Root of f(x) = 0 by the Newton-Raphson method.
%
%   [root, iter, history, flag] = newton_raphson(f, df, x0, es, maxit)
%
%   INPUTS
%     f     : function handle of the equation,   e.g. f  = @(x) x.^2 - 2
%     df    : function handle of the derivative, e.g. df = @(x) 2*x
%     x0    : initial guess of the root
%     es    : stopping criterion, approximate relative error in %  (default 1e-4 %)
%     maxit : maximum number of iterations allowed                 (default 50)
%
%   OUTPUTS
%     root    : estimated root (last iterate)
%     iter    : number of iterations performed
%     history : one row per iteration ->  [ i, x_i, f(x_i), f'(x_i), x_(i+1), ea(%) ]
%     flag    :  1  converged (ea <= es)
%                0  maximum number of iterations reached, not converged
%               -1  derivative became zero (division by zero, method fails)
%               -2  f or f' became complex / NaN / Inf (iterate left the
%                   domain of the function or the method diverged)
%
%   METHOD
%     The tangent to f(x) at the current estimate x_i is extended to the
%     x-axis. The point where it crosses the axis is the new estimate:
%
%            x_(i+1) = x_i - f(x_i) / f'(x_i)
%
%     The approximate relative error is
%
%            ea = | (x_(i+1) - x_i) / x_(i+1) | * 100 %
%
%     and iterations stop when ea <= es.
%
%   CEE 215 Project - Newton-Raphson method

% ---------- default values for optional inputs -------------------------
if nargin < 4 || isempty(es),    es    = 1e-4; end
if nargin < 5 || isempty(maxit), maxit = 50;   end

% ---------- initialisation ---------------------------------------------
xr      = x0;            % current estimate of the root
iter    = 0;             % iteration counter
ea      = 100;           % approximate relative error (%)
flag    = 0;             % assume "not converged" until proven otherwise
history = zeros(0, 6);   % iteration record (grows by one row per step)

% ---------- Newton-Raphson iterations ----------------------------------
while true
    fx  = f(xr);         % function value at the current estimate
    dfx = df(xr);        % slope of the tangent at the current estimate

    % (a) the function must be real and finite at the current estimate
    if ~isreal(fx) || ~isreal(dfx) || ~isfinite(fx) || ~isfinite(dfx)
        flag = -2;
        break;
    end

    % (b) a horizontal tangent never meets the x-axis -> method fails
    if dfx == 0
        flag = -1;
        break;
    end

    % (c) Newton-Raphson formula
    xold = xr;
    xr   = xold - fx / dfx;
    iter = iter + 1;

    % (d) approximate relative error
    if xr ~= 0
        ea = abs((xr - xold) / xr) * 100;
    end

    history(iter, :) = [iter, xold, fx, dfx, xr, ea];

    % (e) stopping checks
    if ~isfinite(xr)
        flag = -2;           % iterate blew up to Inf/NaN (divergence)
        break;
    end
    if ea <= es || fx == 0
        flag = 1;            % converged
        break;
    end
    if iter >= maxit
        flag = 0;            % ran out of iterations
        break;
    end
end

root = xr;
end
