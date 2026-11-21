function show_case(history, flag, root, maxRows)
% SHOW_CASE  Print a Newton-Raphson iteration table and its outcome.
%            Used by NR_Limitations_Demo.m
%
%   show_case(history, flag, root, maxRows)
%     history, flag, root : outputs of newton_raphson.m
%     maxRows             : largest number of rows printed (default 10)
%
%   CEE 215 Project - Newton-Raphson method

if nargin < 4, maxRows = 10; end
nr = size(history, 1);

if nr == 0
    fprintf('  (no iteration could be performed)\n');
else
    fprintf('  Iter         x_i          f(x_i)         f''(x_i)          x_i+1        ea (%%)\n');
    for i = 1:min(nr, maxRows)
        fprintf('  %3d  %13.6g  %14.6g  %14.6g  %13.6g  %11.3g\n', real(history(i, :)));
    end
    if nr > maxRows
        fprintf('  ... (%d more rows not shown)\n', nr - maxRows);
    end
end

switch flag
    case 1
        fprintf('  RESULT : converged to %.8g in %d iterations\n', root, nr);
    case 0
        fprintf('  RESULT : NOT converged after %d iterations (last value %s)\n', nr, num2str(root, 6));
    case -1
        fprintf('  RESULT : FAILED - derivative is zero at x = %s (horizontal tangent)\n', num2str(root, 6));
    case -2
        fprintf('  RESULT : FAILED - f(x) is complex, NaN or Inf at x = %s\n', num2str(root, 6));
end
end
