function [n str] = sigdig(x)
% sigdig calculates the position(s) of the first significant digit(s) of x.
% If x is an array, sigdif acts elementwise. The output is returned in the
% format n: 10^n = place of 1st significant digit of x.

% INPUTS
% x = number, vector, or array of numbers

% OUTPUTS
% n = position of 1st significant digit(s) of (elements of) x; negative n
%   indicates a position after the decimal point, 0 or positive before
% str = string to assist formatting of numbers with sprintf/fprintf to be
%   used e.g. in the format: sprintf(['x = ',str],x)

n = floor(log10(roundsd(x,1)));
if n<0
  str = ['%.',num2str(abs(n)),'f'];
else
  str = ['%',num2str(n),'.f'];
end

end

