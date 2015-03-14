function [indices,subsets] = GenSubsetIndices(data_dim,order)

% Generates two cell arrays of length order. The jth element of subsets is
% a matrix containing binary vectors corresponding to all subsets of weight
% j of a set of size data_dim. That is, each column of the matrix is a
% binary vector of weight j and length data_dim. The jth element of indices
% is also an array. Each column contains the indices of the nonzero
% elements of the corresponding column of subsets.
% 
% Written by: DJ Strouse
% Last updated: May 17, 2013 by DJ Strouse
% Part of: MaxEnt code suite
% Used by: fit_MaxEnt.m
% 
% INPUTS
% data_dim [=] positive integer
% order [=] positive integer
% 
% OUTPUTS
% indices [=] order X 1 cell array
% indices{j} [=] j X (data_dim choose j)
% subsets [=] order X 1 cell array
% subsets{j} [=] data_dim X (data_dim choose j)

subsets = cell(order,1); % each entry [=] data_dim X Nsubsets of weight
indices = cell(order,1); % each entry [=] weight X Nsubsets of weight
num_subsets = 2^data_dim;

% loop over all binary vectors
for i = 1:num_subsets-1
  index_vec = int2bin(i,data_dim)';
  weight = sum(index_vec);
  % add them to appropriate list of subsets
  if weight < order+1
    subsets{weight} = cat(2,subsets{weight},index_vec);
    indices{weight} = cat(2,indices{weight},find(index_vec));
  end
  clear index_vec weight;
end
clear i;

end