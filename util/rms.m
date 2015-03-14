function [rmsx] = rms(X,dim)
% rms calculate the root-mean-square of the columns of input array X
% if X is a (row or column vector), the rms of that vector is returned

if ~exist('dim','var')
  if size(X,2)>1 && size(X,1)==1 % if row vector
    dim = 2; % take rms of row vector
  else
    dim = 1; % otherwise do the columns by default
  end
end

rmsx = sqrt(mean(X.^2,dim));

end

