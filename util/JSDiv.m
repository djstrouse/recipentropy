function dist = JSDiv(P,Q)
% Jensen-Shannon divergence of two probability distributions
% P and Q are automatically normalised to have the sum of one on rows 

if size(P)~=size(Q)
    error('P and Q should be the same size!');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

% normalizing the P and Q
Q = Q ./sum(Q);
P = P ./sum(P);

M = 0.5.*(P + Q);

dist = 0.5 * KLDiv(P,M) + 0.5 * KLDiv(Q,M);