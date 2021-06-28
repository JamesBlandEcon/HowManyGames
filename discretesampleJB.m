function [X] = discretesampleJB(P)
% samples from the discrete distribution P
% each row of P is another distribution
[N,K] = size(P);
u = repmat(rand(N,1),[1,K]);
cP = cumsum(P,2);

LLim  = [zeros(N,1) cP(:,1:(K-1))];
RLim = cP;

Bin = (u>=LLim) & (u<RLim);

NN = repmat(1:K,[N,1]);

X = sum(NN.*Bin,2);




end

