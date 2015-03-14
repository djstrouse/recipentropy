function [Asym] = symmetrize(A)

% Takes an upper triangular matrix and copies the upper triangle
% to the lower triangle, creating a symmetric matrix.
% 
% Written by: DJ Strouse
% Last updated: May 17, 2013 by DJ Strouse
% Part of: MaxEnt code suite
% Used by: fit_MaxEnt.m
%
% INPUTS
% A = upper triangular matrix
%
% OUTPUTS
% Asym = symmetrized A

% check if square
if ~issquare(A)
  error('Non-square matrix cannot be symmetrized!')
end
if sum(sum(A.*A'.*not(eye(size(A,1))),1),2)>0
  error('Transposed elements simultaneously nonzero!')
end

Asym = A+not(eye(size(A,1))).*A';

end

