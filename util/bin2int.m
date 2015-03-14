function [int] = bin2int(binvec)
% Converts a binary vector to its integer equivalent

% confirm binvec is binary vector
if ~isvector(binvec)
  error('Input is not a vector!')
elseif ~isbinary(binvec)
  error('Input is not binary!')
end

% init
data_dim = length(binvec);

% construct binary powers vector
powers = (2*ones(1,data_dim)).^(0:data_dim-1);

% calculate integer value of binvec
int = dot(binvec,powers);

end

