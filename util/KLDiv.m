function dist = KLDiv(P,Q)
% KLDiv calculates the Kullback-Leibler divergence between two discrete
% probability distributions P and Q.
% P and Q  are automatically normalised to have the sum of one on rows.

if size(P)~=size(Q)
  error('P and Q should be the same size.')
end
if sum(~isfinite(P(:))) || sum(~isfinite(Q(:)))
  error('the inputs contain non-finite values.') 
end
if sum(P<0)>0
  error('P contains negative entries.')
end
if sum(Q<0)>0
  error('Q contains negative entries.')
end

% normalizing the P and Q
Q = Q./sum(Q);
P = P./sum(P);

% calculate divergence
R = log(P./Q);
R(P==0) = 0; % ensures 0*log(0)=0
dist = P'*R;