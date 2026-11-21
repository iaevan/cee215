function val = read_input(prompt, default, rule, range)
% READ_INPUT  Ask the user for one number, with a default value and a check.
%
%   val = read_input(prompt, default, rule)
%   val = read_input(prompt, default, 'int', [lo hi])
%
%   prompt  : text shown to the user
%   default : value used when the user just presses Enter
%   rule    : 'pos'    -> value must be > 0
%             'nonneg' -> value must be >= 0
%             'int'    -> value must be a whole number between range(1) and range(2)
%
%   The question is repeated until an acceptable value is typed.
%   Simple expressions such as 1/2000 or 0.5e-3 are also accepted.
%
%   CEE 215 Project - Newton-Raphson method

if nargin < 4, range = [-Inf Inf]; end

while true
    str = input(sprintf('%s [default = %g]: ', prompt, default), 's');

    % Enter only -> take the default value
    if isempty(strtrim(str))
        val = default;
        return;
    end

    val = str2num(str); %#ok<ST2NM>   (str2num so that 1/2000 is accepted)

    if isempty(val) || ~isscalar(val) || ~isreal(val) || ~isfinite(val)
        fprintf('   >> "%s" is not a valid number. Please try again.\n', str);
        continue;
    end

    switch rule
        case 'pos'
            if val > 0, return; end
            fprintf('   >> The value must be greater than zero.\n');
        case 'nonneg'
            if val >= 0, return; end
            fprintf('   >> The value must not be negative.\n');
        case 'int'
            if val == round(val) && val >= range(1) && val <= range(2), return; end
            fprintf('   >> Please enter a whole number from %d to %d.\n', range(1), range(2));
        otherwise
            return;
    end
end
end
