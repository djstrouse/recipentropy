function [fout] = f(data,order,subsets)

% Builds a design matrix for the purpose of fitting a maxent model.
% 
% Written by: DJ Strouse
% Last updated: May 17, 2013 by DJ Strouse
% Part of: MaxEnt code suite
% Used by: fit_MaxEnt.m
% 
% INPUTS
% data [=] data_dim X samp_size
% order [=] positive integer
% subsets [=] order X 1 cell array
% subsets{j} [=] data_dim X (data_dim choose j)
% 
% OUTPUTS
% fout [=] (sum over i=1:order (data_dim choose i)) X samp_size
%          [term in paren above = num param in maxent of order]

data_dim = size(data,1);
samp_size = size(data,2);

if order>data_dim
	error('Order is larger than data set!')
end

if order<1
	error('Order must be a positive integer!')
end

fout_length = data_dim;
for i = 2:order
  fout_length = fout_length+nchoosek(data_dim,i);
end
clear i;
% fout_length counts the number of parameters for order maxent model

% initialize fout with data
fout = zeros(fout_length,samp_size);
fout(1:data_dim,:) = data;

% loop over orders
terms_added = data_dim;
for i = 2:order
  rows_to_add = round(subsets{i}'*data);
  % The j,kth element of rows_to_add is now equal to the sum of terms
  % in the jth subset of size i of the kth sample (i.e. number of terms
  % equal to 1 for binary activity vectors). If that sum is equal to
  % order i (i.e. all terms equal to 1), the product will be 1.
  % Otherwise, it will be zero.
  rows_to_add = rows_to_add==i;
  num_orderi_terms = size(rows_to_add,1);
  fout(terms_added+1:terms_added+num_orderi_terms,:) = rows_to_add;
  terms_added = terms_added+num_orderi_terms;
end
clear i;

end