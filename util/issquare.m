function [answer] = issquare(X)
% issquare tests whether a matrix X is square

% INPUTS
% X = matrix

% OUTPUTS
% answer = boolean indicating whether X is a square matrix

if size(X,1)==size(X,2)
  answer = true;
else
  answer = false;
end

end

