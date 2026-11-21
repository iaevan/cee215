%% ========================================================================
%  CEE 215 : Numerical Analysis and Computer Programming (Section 1)
%  North South University - Department of Civil and Environmental Engg.
%
%  PROJECT : Normal Depth of Uniform Flow in a Trapezoidal Open Channel
%            using the Newton-Raphson Method
%
%  THE PROBLEM
%    For a prismatic trapezoidal channel (rectangular and triangular
%    channels are special cases) carrying a discharge Q in uniform flow,
%    Manning's equation is
%
%         Q = (k/n) * A * R^(2/3) * S0^(1/2)
%
%    A, P and R all depend on the depth y, so the equation cannot be
%    rearranged to give y directly. The script solves
%
%         f(y) = (k/n) * S0^(1/2) * A^(5/3) / P^(2/3) - Q = 0
%
%    for the normal depth yn with the Newton-Raphson method. The critical
%    depth yc (from A^3/T = Q^2/g) is also found with Newton-Raphson, and
%    the velocity, Froude number, flow regime and slope type are reported.
%
%  FILES
%    NR_Channel_Main.m      <- RUN THIS FILE
%    newton_raphson.m       general Newton-Raphson solver (function)
%    read_input.m           keyboard input with a default value and checks
%    NR_Limitations_Demo.m  examples showing where the method/script fails
%
%  OUTPUT
%    Iteration tables and results in the Command Window, three figures,
%    and a text summary, all saved in the "results" folder.
%% ========================================================================

clc; clear; close all;

%% ------------------------------------------------------------------------
%  1. CHOOSE A PROBLEM AND READ THE DATA
%  ------------------------------------------------------------------------
fprintf('=====================================================================\n');
fprintf('   NORMAL DEPTH IN A TRAPEZOIDAL CHANNEL - NEWTON-RAPHSON METHOD\n');
fprintf('   CEE 215 : Numerical Analysis and Computer Programming\n');
fprintf('=====================================================================\n\n');
fprintf('  1 = Example 1 : earthen irrigation canal (SI units)\n');
fprintf('  2 = Example 2 : textbook check, Chow (1959) Example 6.2 (US units)\n');
fprintf('  3 = Enter your own channel data\n\n');
option = read_input('Select an option', 1, 'int', [1 3]);

switch option
    case 1      % trapezoidal earthen irrigation canal
        unitSys = 1;  Q = 25;   b = 6;   z = 1.5;  n = 0.025;  S0 = 0.0004;
    case 2      % Chow, Open-Channel Hydraulics (1959), Example 6.2: yn = 3.36 ft
        unitSys = 2;  Q = 400;  b = 20;  z = 2;    n = 0.025;  S0 = 0.0016;
    case 3      % user's own channel
        fprintf('\n  Unit system:  1 = SI (m, m^3/s)    2 = US customary (ft, ft^3/s)\n');
        unitSys = read_input('Unit system', 1, 'int', [1 2]);
        fprintf('\n  Channel data  (rectangular channel: z = 0,  triangular channel: b = 0)\n');
        Q  = read_input('  Discharge, Q', 25, 'pos');
        b  = read_input('  Bottom width, b', 6, 'nonneg');
        z  = read_input('  Side slope, z (z horizontal : 1 vertical)', 1.5, 'nonneg');
        while b == 0 && z == 0
            fprintf('   >> b and z cannot both be zero (the channel would have no width).\n');
            b = read_input('  Bottom width, b', 6, 'nonneg');
            z = read_input('  Side slope, z (z horizontal : 1 vertical)', 1.5, 'nonneg');
        end
        n  = read_input('  Manning roughness coefficient, n', 0.025, 'pos');
        S0 = read_input('  Longitudinal bed slope, S0', 0.0004, 'pos');
end

% ---- unit system constants ----
if unitSys == 1
    k = 1.0;    g = 9.81;  L = 'm';   Qu = 'm^3/s';   unitName = 'SI';
else
    k = 1.486;  g = 32.2;  L = 'ft';  Qu = 'ft^3/s';  unitName = 'US customary';
end

% ---- plausibility warnings (the program still runs) ----
if n < 0.008 || n > 0.20
    fprintf('\n  WARNING: n = %g is outside the usual range 0.008 - 0.20. Please check.\n', n);
end
if S0 > 0.1
    fprintf('\n  WARNING: S0 = %g is very steep for Manning''s equation. Please check.\n', S0);
end

% ---- name of the cross-section ----
if z == 0
    shape = 'Rectangular';
elseif b == 0
    shape = 'Triangular';
else
    shape = 'Trapezoidal';
end

%% ------------------------------------------------------------------------
%  2. CHANNEL GEOMETRY AND THE EQUATIONS TO BE SOLVED
%  ------------------------------------------------------------------------
m  = sqrt(1 + z^2);              % length of one sloping side per unit depth
C  = k * sqrt(S0) / n;           % constant part of Manning's equation

A  = @(y) (b + z*y) .* y;        % flow area,             A = (b + z y) y
P  = @(y) b + 2*m*y;             % wetted perimeter,      P = b + 2 y sqrt(1+z^2)
T  = @(y) b + 2*z*y;             % top width,             T = b + 2 z y  (= dA/dy)

% Normal depth:   f(y) = C A^(5/3) P^(-2/3) - Q = 0
f  = @(y) C * A(y).^(5/3) ./ P(y).^(2/3) - Q;

% Derivative by the product and chain rules, using dA/dy = T and dP/dy = 2m:
%   f'(y) = C [ (5/3) A^(2/3) T P^(-2/3)  -  (2/3) A^(5/3) (2m) P^(-5/3) ]
df = @(y) C * ( (5/3) * A(y).^(2/3) .* T(y) ./ P(y).^(2/3) ...
              - (2/3) * A(y).^(5/3) .* (2*m) ./ P(y).^(5/3) );

% Critical depth (Froude number = 1):   fc(y) = A^3/T - Q^2/g = 0
%   fc'(y) = 3 A^2 - 2 z A^3 / T^2
fc  = @(y) A(y).^3 ./ T(y) - Q^2/g;
dfc = @(y) 3*A(y).^2 - 2*z*A(y).^3 ./ T(y).^2;

%% ------------------------------------------------------------------------
%  3. NEWTON-RAPHSON SETTINGS
%  ------------------------------------------------------------------------
% Automatic first guess from a simple shape whose depth has a closed form
if b > 0
    y0auto = (Q*n / (k*b*sqrt(S0)))^(3/5);              % wide rectangle: A = b y, R = y
else
    y0auto = (Q / (C * z^(5/3) / (2*m)^(2/3)))^(3/8);    % triangle: A = z y^2, P = 2 m y
end

fprintf('\n  Initial guess:  1 = automatic (wide-channel approximation)   2 = type it\n');
guessType = read_input('Initial guess option', 1, 'int', [1 2]);
if guessType == 1
    y0 = y0auto;
    guessText = 'automatic';
    fprintf('  Automatic initial guess y0 = %.4f %s\n', y0, L);
else
    y0 = read_input(sprintf('  Initial guess, y0 (%s)', L), round(y0auto*100)/100, 'pos');
    guessText = 'user input';
end

nsig  = read_input('  Number of significant figures required', 6, 'int', [1 12]);
maxit = read_input('  Maximum number of iterations', 50, 'int', [1 1000]);
es    = 0.5 * 10^(2 - nsig);     % Scarborough (1966): result correct to nsig figures
fprintf('  Stopping criterion: es = 0.5 x 10^(2-%d) = %g %%\n', nsig, es);

%% ------------------------------------------------------------------------
%  4. SOLVE FOR THE NORMAL DEPTH  (Newton-Raphson)
%  ------------------------------------------------------------------------
[yn, itN, histN, flagN] = newton_raphson(f, df, y0, es, maxit);

if flagN ~= 1
    fprintf('\n  Iterations performed before stopping:\n');
    fprintf('  Iter        y_i            f(y_i)         f''(y_i)          y_i+1\n');
    for i = 1:itN
        fprintf('  %3d  %14.6g  %14.6g  %14.6g  %14.6g\n', histN(i, 1:5));
    end
    fprintf('\n  *** The Newton-Raphson method did NOT find the normal depth. ***\n');
    switch flagN
        case 0
            fprintf('  Reason: no convergence within %d iterations.\n', maxit);
        case -1
            fprintf('  Reason: f''(y) = 0 at y = %g, the tangent is horizontal.\n', yn);
        case -2
            fprintf('  Reason: f(y) is not real/finite at y = %g (depth left the physical range).\n', yn);
    end
    fprintf('  Try another initial guess, more iterations, or check the data.\n');
    return;
end

%% ------------------------------------------------------------------------
%  5. SOLVE FOR THE CRITICAL DEPTH  (Newton-Raphson again)
%  ------------------------------------------------------------------------
if b > 0
    yc0 = ((Q/b)^2 / g)^(1/3);          % critical depth of a rectangle of width b
else
    yc0 = (2*Q^2 / (g*z^2))^(1/5);      % exact critical depth of a triangle
end
[yc, itC, histC, flagC] = newton_raphson(fc, dfc, yc0, es, maxit);
if flagC ~= 1
    fprintf('\n  NOTE: critical depth did not converge (flag %d); it is not reported.\n', flagC);
    yc = NaN;
end

%% ------------------------------------------------------------------------
%  6. HYDRAULIC PROPERTIES AT NORMAL DEPTH
%  ------------------------------------------------------------------------
An = A(yn);   Pn = P(yn);   Tn = T(yn);
Rn = An / Pn;                    % hydraulic radius
V  = Q / An;                     % mean velocity
D  = An / Tn;                    % hydraulic depth
Fr = V / sqrt(g*D);              % Froude number
E  = yn + V^2 / (2*g);           % specific energy

if abs(Fr - 1) < 0.01
    regime = 'Critical flow (Fr = 1)';
elseif Fr < 1
    regime = 'Subcritical / tranquil flow (Fr < 1)';
else
    regime = 'Supercritical / rapid flow (Fr > 1)';
end

if isnan(yc)
    slopeType = 'not available';
elseif abs(yn - yc) / yc < 0.01
    slopeType = 'Critical slope, C  (yn = yc)';
elseif yn > yc
    slopeType = 'Mild slope, M  (yn > yc)';
else
    slopeType = 'Steep slope, S  (yn < yc)';
end

%% ------------------------------------------------------------------------
%  7. VERIFICATION OF THE ANSWER
%  ------------------------------------------------------------------------
Qcheck = C * An * Rn^(2/3);                      % put yn back into Manning
errQ   = abs(Qcheck - Q) / Q * 100;
yfz    = fzero(f, [1e-6*yn, 10*yn]);             % MATLAB built-in solver
errFz  = abs(yn - yfz) / yfz * 100;

%% ------------------------------------------------------------------------
%  8. PRINT RESULTS  (to the Command Window and to a text file)
%  ------------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
resDir = fullfile(scriptDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

fid  = fopen(fullfile(resDir, 'NR_Channel_Results.txt'), 'w');
outs = 1;                          % 1 = Command Window
if fid > 0, outs = [1 fid]; end    % also write to the file if it opened

for out = outs
    fprintf(out, '\n=====================================================================\n');
    fprintf(out, ' NORMAL DEPTH BY THE NEWTON-RAPHSON METHOD - RESULTS\n');
    fprintf(out, '=====================================================================\n');
    fprintf(out, ' INPUT DATA (%s units)\n', unitName);
    fprintf(out, '   Channel shape ......................... %s\n', shape);
    fprintf(out, '   Discharge, Q .......................... %.4f %s\n', Q, Qu);
    fprintf(out, '   Bottom width, b ....................... %.4f %s\n', b, L);
    fprintf(out, '   Side slope, z (H:V) ................... %.4f\n', z);
    fprintf(out, '   Manning''s n ........................... %.4f\n', n);
    fprintf(out, '   Bed slope, S0 ......................... %.6f\n', S0);
    fprintf(out, ' SOLVER SETTINGS\n');
    fprintf(out, '   Initial guess, y0 ..................... %.4f %s (%s)\n', y0, L, guessText);
    fprintf(out, '   Significant figures ................... %d  (es = %g %%)\n', nsig, es);
    fprintf(out, '   Maximum iterations .................... %d\n', maxit);

    fprintf(out, '\n NEWTON-RAPHSON ITERATIONS FOR NORMAL DEPTH, yn\n');
    fprintf(out, ' Iter     y_i (%s)         f(y_i)        f''(y_i)     y_i+1 (%s)       ea (%%)\n', L, L);
    for i = 1:itN
        fprintf(out, ' %3d  %13.8f  %13.5e  %13.5e  %13.8f  %12.4e\n', histN(i, :));
    end

    if ~isnan(yc)
        fprintf(out, '\n NEWTON-RAPHSON ITERATIONS FOR CRITICAL DEPTH, yc\n');
        fprintf(out, ' Iter     y_i (%s)        fc(y_i)       fc''(y_i)     y_i+1 (%s)       ea (%%)\n', L, L);
        for i = 1:itC
            fprintf(out, ' %3d  %13.8f  %13.5e  %13.5e  %13.8f  %12.4e\n', histC(i, :));
        end
    end

    fprintf(out, '\n HYDRAULIC RESULTS\n');
    fprintf(out, '   Normal depth, yn ...................... %.6f %s  (%d iterations)\n', yn, L, itN);
    if ~isnan(yc)
        fprintf(out, '   Critical depth, yc .................... %.6f %s  (%d iterations)\n', yc, L, itC);
    end
    fprintf(out, '   Flow area, A .......................... %.4f %s^2\n', An, L);
    fprintf(out, '   Wetted perimeter, P ................... %.4f %s\n', Pn, L);
    fprintf(out, '   Hydraulic radius, R = A/P ............. %.4f %s\n', Rn, L);
    fprintf(out, '   Top width, T .......................... %.4f %s\n', Tn, L);
    fprintf(out, '   Mean velocity, V = Q/A ................ %.4f %s/s\n', V, L);
    fprintf(out, '   Hydraulic depth, D = A/T .............. %.4f %s\n', D, L);
    fprintf(out, '   Froude number, Fr ..................... %.4f\n', Fr);
    fprintf(out, '   Specific energy, E .................... %.4f %s\n', E, L);
    fprintf(out, '   Flow regime ........................... %s\n', regime);
    fprintf(out, '   Channel slope type .................... %s\n', slopeType);

    fprintf(out, '\n VERIFICATION\n');
    fprintf(out, '   Q recomputed from yn .................. %.6f %s  (error %.2e %%)\n', Qcheck, Qu, errQ);
    fprintf(out, '   MATLAB fzero solution ................. %.6f %s  (difference %.2e %%)\n', yfz, L, errFz);
    fprintf(out, '=====================================================================\n');
end
if fid > 0, fclose(fid); end

%% ------------------------------------------------------------------------
%  9. FIGURES
%  ------------------------------------------------------------------------
% ---- Figure 1: f(y) with the Newton-Raphson tangent steps ----
fig1 = figure('Name', 'Newton-Raphson on f(y)', 'NumberTitle', 'off');
yAll = [histN(:, 2); histN(:, 5)];
yMax = max(1.6*yn, 1.1*max(yAll));
yy   = linspace(0, yMax, 400);
hF   = plot(yy, f(yy), 'b-', 'LineWidth', 2);  hold on;  grid on;
plot([0 yMax], [0 0], 'k-', 'LineWidth', 0.8);
nShow = min(itN, 6);                                % first few steps are enough
for i = 1:nShow
    yi = histN(i, 2);   fi = histN(i, 3);   yi1 = histN(i, 5);
    hT = plot([yi yi1], [fi 0], 'r--', 'LineWidth', 1.3);         % tangent line
    hV = plot([yi1 yi1], [0 f(yi1)], 'k:', 'LineWidth', 1.2);     % back to curve
    hP = plot(yi, fi, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    % label the point only if it is not on top of the previous label
    if i == 1 || abs(yi - histN(i-1, 2)) > 0.04*yMax
        text(yi, fi, sprintf('  y_%d', i-1), 'FontSize', 11, 'VerticalAlignment', 'top');
    end
end
hR = plot(yn, 0, 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
xlabel(sprintf('Depth, y (%s)', L));
ylabel(sprintf('f(y)  (%s)', Qu));
title(sprintf('Newton-Raphson: y_n = %.4f %s after %d iterations', yn, L, itN));
legend([hF hT hV hP hR], {'f(y) = Q_{Manning}(y) - Q', 'Tangent at y_i', ...
       'Back to the curve', 'Iterates y_i', 'Root, y_n'}, 'Location', 'northwest');

% ---- Figure 2: convergence history ----
fig2 = figure('Name', 'Convergence history', 'NumberTitle', 'off');
subplot(1, 2, 1);
semilogy(histN(:, 1), histN(:, 6), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;  grid on;
semilogy([0.5 itN+0.5], [es es], 'r--', 'LineWidth', 1.2);
xlim([0.5 itN+0.5]);
xlabel('Iteration, i');  ylabel('\epsilon_a  (%)');
title('Approximate relative error');
legend('\epsilon_a', '\epsilon_s (stopping criterion)', 'Location', 'northeast');

subplot(1, 2, 2);
resid = max(abs(f(histN(:, 5))), 1e-16);     % floor so that log scale works
semilogy(histN(:, 1), resid, 'ms-', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
grid on;  xlim([0.5 itN+0.5]);
xlabel('Iteration, i');  ylabel(sprintf('|f(y_{i+1})|  (%s)', Qu));
title('Residual of Manning''s equation');

% ---- Figure 3: channel cross-section ----
fig3 = figure('Name', 'Channel cross-section', 'NumberTitle', 'off');
H  = 1.3 * max(yn, yc);                      % bank height drawn (max ignores NaN)
xL = @(y) -(b/2 + z*y);                      % left bank x-coordinate at depth y
xR = @(y)   b/2 + z*y;                       % right bank x-coordinate at depth y
hW = fill([xL(yn) -b/2 b/2 xR(yn)], [yn 0 0 yn], [0.60 0.80 0.97], 'EdgeColor', 'none');
hold on;
hC = plot([xL(H) -b/2 b/2 xR(H)], [H 0 0 H], 'k-', 'LineWidth', 2.5);
hN = plot([xL(yn) xR(yn)], [yn yn], 'b-', 'LineWidth', 2);
legItems = [hW hC hN];
legText  = {'Water at uniform flow', 'Channel boundary', sprintf('Normal depth y_n = %.3f %s', yn, L)};
if ~isnan(yc)
    hK = plot([xL(yc) xR(yc)], [yc yc], 'r--', 'LineWidth', 2);
    legItems(end+1) = hK;
    legText{end+1}  = sprintf('Critical depth y_c = %.3f %s', yc, L);
end
axis equal;  grid on;
xlim([1.08*xL(H) 1.08*xR(H)]);          % whole section, true scale (axis equal)
ylim([-0.08*H 2.3*H]);                 % head-room for the legend
xlabel(sprintf('Width (%s)', L));  ylabel(sprintf('Height above bed (%s)', L));
title(sprintf('%s channel: Q = %g %s, V = %.2f %s/s, Fr = %.2f', shape, Q, Qu, V, L, Fr));
legend(legItems, legText, 'Location', 'north');

% ---- save the figures ----
print(fig1, '-dpng', '-r150', fullfile(resDir, 'Fig1_NR_iterations.png'));
print(fig2, '-dpng', '-r150', fullfile(resDir, 'Fig2_convergence.png'));
print(fig3, '-dpng', '-r150', fullfile(resDir, 'Fig3_cross_section.png'));

fprintf('\n Figures and NR_Channel_Results.txt saved in:\n   %s\n\n', resDir);
