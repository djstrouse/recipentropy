function [negL,delta] = NegLogLik(data_avgd_f,num_samples,data_dim,order,subsets,lambda)

% MaxEnt log likelihood of data using parameters lambda
% assumes data is the output of f applied to a data matrix

num_patterns = 2^data_dim;

% calculate Z and dist_avgd_f
Z = 0;
dist_avgd_f = 0;
for i = 1:num_patterns % can this be written as a matrix multiplication?
  pattern = int2bin(i-1,data_dim);
  f_pattern = f(pattern',order,subsets);
  p = exp(lambda'*f_pattern);
  Z = Z+p;
  dist_avgd_f = dist_avgd_f+p*f_pattern';
  clear p;
end
clear i;

% normalize
dist_avgd_f = dist_avgd_f/Z;

% calculate neg log lik

L1 = -num_samples*log(Z);
L2 = num_samples*lambda'*data_avgd_f;
negL = -L1-L2;

% calculate the gradient

delta1 = num_samples*data_avgd_f';
delta2 = -num_samples*dist_avgd_f;

% flip the sign since you are using the NEGATIVE log lik
delta = -delta1-delta2;

end

