%% ========================================================================
%  NR_Limitations_Demo.m
%  CEE 215 Project - LIMITATIONS of the Newton-Raphson method and of the
%  normal-depth script, each shown with a worked example.
%
%  Cases
%    1. Effect of the initial guess on the channel problem
%    2. Initial guess y0 = 0 (dry channel)        -> f'(y) = 0, method fails
%    3. Negative initial guess                    -> complex numbers, fails
%    4. Horizontal channel, S0 = 0                -> no normal depth exists
%    5. Oscillation      f(x) = x^3 - 2x + 2,          x0 = 0
%    6. Divergence       f(x) = atan(x),               x0 = 1.5
%    7. Multiple root    f(x) = (x-1)^2 (x-3)          -> slow convergence
%    8. Circular culvert flowing nearly full      -> two depths give the
%       same Q; the answer depends on the initial guess
%
%  Needs: newton_raphson.m, show_case.m
%% ========================================================================

clc; clear; close all;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
resDir = fullfile(scriptDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

nextCase = @() input('\n  Press Enter for the next case ... ', 's');

es    = 5e-5;      % 6 significant figures (Scarborough criterion, %)
maxit = 50;

% ---- Example 1 channel from NR_Channel_Main.m ----
Q = 25;  b = 6;  z = 1.5;  n = 0.025;  S0 = 0.0004;  k = 1.0;
m  = sqrt(1 + z^2);
C  = k * sqrt(S0) / n;
A  = @(y) (b + z*y) .* y;
P  = @(y) b + 2*m*y;
T  = @(y) b + 2*z*y;
f  = @(y) C * A(y).^(5/3) ./ P(y).^(2/3) - Q;
df = @(y) C * ( (5/3) * A(y).^(2/3) .* T(y) ./ P(y).^(2/3) ...
              - (2/3) * A(y).^(5/3) .* (2*m) ./ P(y).^(5/3) );

fprintf('=====================================================================\n');
fprintf('   LIMITATIONS OF THE NEWTON-RAPHSON METHOD AND OF THE SCRIPT\n');
fprintf('=====================================================================\n');
fprintf('  Channel used in Cases 1-4: Q = %g m^3/s, b = %g m, z = %g, n = %g, S0 = %g\n', ...
        Q, b, z, n, S0);
fprintf('  Correct normal depth (from NR_Channel_Main.m): yn = 2.403895 m\n');

%% ------------------------------------------------------------------------
%  CASE 1 : effect of the initial guess
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 1 : Effect of the initial guess y0 on the number of iterations\n');
fprintf('---------------------------------------------------------------------\n');
y0List = [0.01 0.1 0.5 1 2 2.4 5 10 50 100];
nG     = numel(y0List);
itList = zeros(1, nG);
y1List = zeros(1, nG);
fprintf('     y0 (m)      first step y1 (m)   iterations      root (m)\n');
for j = 1:nG
    [r, it, h, fl] = newton_raphson(f, df, y0List(j), es, maxit);
    itList(j) = it;
    y1List(j) = h(1, 5);
    fprintf('   %8.2f      %14.4f       %5d         %10.6f\n', y0List(j), y1List(j), it, r);
end
fprintf('\n  LESSON: every positive guess reaches yn here, but a poor guess costs\n');
fprintf('  many more iterations. A guess near zero gives a nearly flat tangent,\n');
fprintf('  so the first step overshoots to a very large depth.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 2 : y0 = 0  -> f'(0) = 0
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 2 : Initial guess y0 = 0 (dry channel)\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('  f(0) = %g,  f''(0) = %g\n', f(0), df(0));
[r, it, h, fl] = newton_raphson(f, df, 0, es, maxit);
show_case(h, fl, r);
fprintf('  LESSON: at y = 0 the area is zero, so f''(0) = 0. The tangent is\n');
fprintf('  horizontal and never meets the axis: y1 = y0 - f/f'' needs a division\n');
fprintf('  by zero. The main script only accepts y0 > 0 for this reason.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 3 : negative initial guess -> complex numbers
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 3 : Negative initial guess y0 = -1 m\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('  A(-1) = %g m^2  ->  A^(5/3) = %s  (complex)\n', A(-1), num2str(A(-1)^(5/3), 5));
[r, it, h, fl] = newton_raphson(f, df, -1, es, maxit);
show_case(h, fl, r);
fprintf('  LESSON: a negative depth has no physical meaning; the fractional\n');
fprintf('  powers in Manning''s equation turn complex and the method stops.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 4 : horizontal bed, S0 = 0
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 4 : Horizontal channel bed, S0 = 0\n');
fprintf('---------------------------------------------------------------------\n');
C0  = k * sqrt(0) / n;                                       % = 0
f0  = @(y) C0 * A(y).^(5/3) ./ P(y).^(2/3) - Q;              % = -Q for every y
df0 = @(y) C0 * ( (5/3) * A(y).^(2/3) .* T(y) ./ P(y).^(2/3) ...
                - (2/3) * A(y).^(5/3) .* (2*m) ./ P(y).^(5/3) );
fprintf('  f(y) = %g for every y,  f''(y) = %g\n', f0(2), df0(2));
[r, it, h, fl] = newton_raphson(f0, df0, 2.4, es, maxit);
show_case(h, fl, r);
fprintf('  LESSON: Manning''s equation gives no uniform flow on a horizontal (or\n');
fprintf('  adverse) slope, so no root exists. The main script therefore accepts\n');
fprintf('  only S0 > 0.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 5 : oscillation
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 5 : Oscillation,  f(x) = x^3 - 2x + 2,  x0 = 0\n');
fprintf('---------------------------------------------------------------------\n');
f5  = @(x) x.^3 - 2*x + 2;
df5 = @(x) 3*x.^2 - 2;
[r5, it5, h5, fl5] = newton_raphson(f5, df5, 0, es, maxit);
show_case(h5, fl5, r5, 6);
[r5b, it5b, h5b, fl5b] = newton_raphson(f5, df5, -2, es, maxit);
fprintf('  Same function from x0 = -2:\n');
show_case(h5b, fl5b, r5b);
fprintf('  LESSON: the iterates jump 0 -> 1 -> 0 -> 1 ... forever. The method\n');
fprintf('  is trapped in a cycle and only the iteration limit stops it.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 6 : divergence
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 6 : Divergence,  f(x) = atan(x),  x0 = 1.5  (true root x = 0)\n');
fprintf('---------------------------------------------------------------------\n');
f6  = @(x) atan(x);
df6 = @(x) 1 ./ (1 + x.^2);
[r6, it6, h6, fl6] = newton_raphson(f6, df6, 1.5, es, maxit);
show_case(h6, fl6, r6, 8);
[r6b, it6b, h6b, fl6b] = newton_raphson(f6, df6, 1.0, es, maxit);
fprintf('  Same function from x0 = 1.0: root x = %g found in %d iterations.\n', r6b, it6b);
fprintf('  LESSON: far from the root the tangent is almost flat, every step\n');
fprintf('  overshoots further and the iterates run away to infinity (until the\n');
fprintf('  slope underflows to zero). For atan(x), |x0| must be below about 1.39.\n');
nextCase();

%% ------------------------------------------------------------------------
%  CASE 7 : multiple (repeated) root
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 7 : Multiple root,  f(x) = (x-1)^2 (x-3) = x^3 - 5x^2 + 7x - 3\n');
fprintf('---------------------------------------------------------------------\n');
f7  = @(x) x.^3 - 5*x.^2 + 7*x - 3;
df7 = @(x) 3*x.^2 - 10*x + 7;
[r7, it7, h7, fl7] = newton_raphson(f7, df7, 0, es, maxit);
fprintf('  Double root x = 1, starting at x0 = 0:\n');
show_case(h7, fl7, r7, 6);
[r7b, it7b, h7b, fl7b] = newton_raphson(f7, df7, 4, es, maxit);
fprintf('  Simple root x = 3, starting at x0 = 4:\n');
show_case(h7b, fl7b, r7b);
fprintf('  LESSON: at a double root f''(x) -> 0 as well, so convergence drops\n');
fprintf('  from quadratic to linear (%d iterations instead of %d).\n', it7, it7b);
nextCase();

%% ------------------------------------------------------------------------
%  CASE 8 : circular culvert - two depths carry the same discharge
%  ------------------------------------------------------------------------
fprintf('\n---------------------------------------------------------------------\n');
fprintf(' CASE 8 : Circular culvert (outside the scope of the main script)\n');
fprintf('---------------------------------------------------------------------\n');
Dc = 1.0;  nc = 0.013;  Sc = 0.002;  Qc = 1.10;
th  = @(y) 2 * acos(1 - 2*y/Dc);                    % central angle (rad)
Acf = @(y) Dc^2/8 * (th(y) - sin(th(y)));           % flow area
Pcf = @(y) Dc * th(y) / 2;                          % wetted perimeter
Qm  = @(y) (1/nc) * sqrt(Sc) * Acf(y).^(5/3) ./ Pcf(y).^(2/3);
f8  = @(y) Qm(y) - Qc;
hd  = 1e-6;                                         % step for the derivative
df8 = @(y) (f8(y + hd) - f8(y - hd)) / (2*hd);      % central difference
fprintf('  D = %g m, n = %g, S0 = %g, Q = %g m^3/s\n', Dc, nc, Sc, Qc);
fprintf('  Discharge when flowing just full, Q(D) = %.4f m^3/s\n', Qm(Dc));
fprintf('  Largest possible discharge Q(0.938 D) = %.4f m^3/s\n', Qm(0.938*Dc));
[r8a, it8a, h8a, fl8a] = newton_raphson(f8, df8, 0.5, es, maxit);
fprintf('\n  Initial guess y0 = 0.50 m:\n');   show_case(h8a, fl8a, r8a);
[r8b, it8b, h8b, fl8b] = newton_raphson(f8, df8, 0.999, es, maxit);
fprintf('\n  Initial guess y0 = 0.999 m:\n');  show_case(h8b, fl8b, r8b);
[r8c, it8c, h8c, fl8c] = newton_raphson(f8, df8, 0.97, es, maxit);
fprintf('\n  Initial guess y0 = 0.97 m:\n');   show_case(h8c, fl8c, r8c);
fprintf('  LESSON: near the crown Q(y) rises and then falls, so two depths\n');
fprintf('  carry the same Q. Newton-Raphson returns whichever root is nearest\n');
fprintf('  the guess and cannot tell the user that a second root exists. A\n');
fprintf('  guess very close to D can also jump above the crown (y > D), where\n');
fprintf('  acos() turns complex. The main script handles only trapezoidal,\n');
fprintf('  rectangular and triangular channels, where Q(y) always increases.\n');

%% ------------------------------------------------------------------------
%  FIGURES
%  ------------------------------------------------------------------------
% ---- Figure L1: iterations versus initial guess ----
figL1 = figure('Name', 'Effect of initial guess', 'NumberTitle', 'off');
semilogx(y0List, itList, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');  hold on;  grid on;
plot([2.403895 2.403895], [0 max(itList)+2], 'r--', 'LineWidth', 1.2);
xlabel('Initial guess, y_0 (m)  [log scale]');  ylabel('Iterations to converge');
title('Case 1: iterations needed versus initial guess');
legend('Newton-Raphson', 'True normal depth y_n = 2.404 m', 'Location', 'north');
ylim([0 max(itList)+2]);

% ---- Figure L2: four failure modes ----
figL2 = figure('Name', 'Failure modes', 'NumberTitle', 'off', 'Position', [100 100 1000 760]);

subplot(2, 2, 1);                                   % (a) oscillation
xx = linspace(-2.2, 1.6, 300);
plot(xx, f5(xx), 'b-', 'LineWidth', 1.8);  hold on;  grid on;
plot([-2.2 1.6], [0 0], 'k-');
for i = 1:4
    plot([h5(i, 2) h5(i, 5)], [h5(i, 3) 0], 'r--', 'LineWidth', 1.2);
    plot([h5(i, 5) h5(i, 5)], [0 f5(h5(i, 5))], 'k:', 'LineWidth', 1.2);
end
plot([0 1], [f5(0) f5(1)], 'ro', 'MarkerFaceColor', 'r');
plot(r5b, 0, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
title('(a) Oscillation: x^3 - 2x + 2, x_0 = 0');  xlabel('x');  ylabel('f(x)');

subplot(2, 2, 2);                                   % (b) divergence
xx = linspace(-12, 12, 400);
plot(xx, f6(xx), 'b-', 'LineWidth', 1.8);  hold on;  grid on;
plot([-12 12], [0 0], 'k-');
for i = 1:min(3, size(h6, 1))
    plot([h6(i, 2) h6(i, 5)], [h6(i, 3) 0], 'r--', 'LineWidth', 1.2);
    plot([h6(i, 5) h6(i, 5)], [0 f6(h6(i, 5))], 'k:', 'LineWidth', 1.2);
    plot(h6(i, 2), h6(i, 3), 'ro', 'MarkerFaceColor', 'r');
end
xlim([-12 12]);
title('(b) Divergence: atan(x), x_0 = 1.5');  xlabel('x');  ylabel('f(x)');

subplot(2, 2, 3);                                   % (c) multiple root
semilogy(h7(:, 1), h7(:, 6), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');  hold on;  grid on;
semilogy(h7b(:, 1), h7b(:, 6), 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
xlabel('Iteration');  ylabel('\epsilon_a (%)');
title('(c) Double root (linear) vs simple root (quadratic)');
legend(sprintf('Double root x = 1 (%d it.)', it7), sprintf('Simple root x = 3 (%d it.)', it7b), ...
       'Location', 'southwest');

subplot(2, 2, 4);                                   % (d) circular culvert
yy = linspace(0.02, Dc, 300);
plot(yy, Qm(yy), 'b-', 'LineWidth', 1.8);  hold on;  grid on;
plot([0 Dc], [Qc Qc], 'k--', 'LineWidth', 1.2);
plot([r8a r8b], [Qc Qc], 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot([0.5 0.999], Qm([0.5 0.999]), 'ro', 'MarkerFaceColor', 'r');
xlabel('Depth, y (m)');  ylabel('Q (m^3/s)');
title('(d) Circular culvert: two depths for one Q');
legend('Manning Q(y)', sprintf('Q = %.2f m^3/s', Qc), 'Two roots found', 'Initial guesses', ...
       'Location', 'northwest');

print(figL1, '-dpng', '-r150', fullfile(resDir, 'FigL1_initial_guess.png'));
print(figL2, '-dpng', '-r150', fullfile(resDir, 'FigL2_failure_modes.png'));

%% ------------------------------------------------------------------------
%  SUMMARY
%  ------------------------------------------------------------------------
fprintf('\n=====================================================================\n');
fprintf(' SUMMARY OF LIMITATIONS\n');
fprintf('=====================================================================\n');
fprintf('  1. Needs a good initial guess; a poor one costs iterations.\n');
fprintf('  2. Fails when f''(x) = 0 (dry channel y0 = 0, horizontal bed).\n');
fprintf('  3. Fails when an iterate leaves the physical domain (y < 0, y > D).\n');
fprintf('  4. Can oscillate or diverge for some functions/guesses.\n');
fprintf('  5. Only linear convergence at multiple roots.\n');
fprintf('  6. Finds one root only - the one nearest the guess.\n');
fprintf('  7. Needs the derivative f''(x), which must be derived by hand.\n');
fprintf('  Figures saved in: %s\n\n', resDir);
