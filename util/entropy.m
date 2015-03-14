function [ent] = entropy(P)
% Takes a matrix whose cols are probability vectors and returns a col
% vector of entropies (in bits)
% P [=] num_probs X num_dist

tolerance = 10^-4;

num_dist = size(P,2);
num_probs = size(P,1);

if max(sum(P,1)<1-tolerance)||max(sum(P,1)>1+tolerance)
  error('At least one probability vector is not normalized!')
end

if isnan(max(sum(P,1)))||isinf(max(sum(P,1)))
  error('At least one probability vector is unnormalizable!')
end

% normalize each col
for i = 1:num_dist
  P(:,i) = P(:,i)./sum(P(:,i));
end
clear i;

% calculate entropy
M = log2(P);
M(isnan(M)) = 0;
M(isinf(M)) = 0;

ent = zeros(num_dist,1);
for i = 1:num_dist
  ent(i,1) = -dot(P(:,i),M(:,i));
end
clear i;

end