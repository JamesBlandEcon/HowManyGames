function [CP,P,CPb,CPn] = RDChoiceProbs(rRD,broadPRi,g,l,ACTCOMB,ACTCOMBPm,D1,D2)

[N,T] = size(ACTCOMB);
u = @(x) x.^(repmat(rRD,[1 T]));

a = 10; c = 56; x = 100; y = 160;

g = repmat(g,[1 1 4]);

P = ones(N,T,4)./4;

num = zeros(size(g));
den = zeros(size(g));


for tt = 2:T
    num = ACTCOMBPm(:,tt-1,:) + g.*num;
    den = 1 + g.*den;
    P(:,tt,:) = num./den;
end

PA = P(:,:,1) + P(:,:,2);
PC = P(:,:,1) + P(:,:,3);

EU = zeros(N,T,4);

EU(:,:,1) = u(a+c+D1+D2); % AC
EU(:,:,2) = PC.*u(a+y+D1+D2)+(1-PC).*u(a+D1+D2); % AD
EU(:,:,3) = PA.*u(x+c+D1+D2)+(1-PA).*u(c+D1+D2); % BC
EU(:,:,4) = P(:,:,1).*u(x+y+D1+D2)+P(:,:,2).*u(x+D1+D2)+P(:,:,3).*u(y+D1+D2)+P(:,:,4).*u(D1+D2);

EUN = zeros(N,T,4);
EUN(:,:,1) = u(a+D1)+u(a+D2);
EUN(:,:,2) = u(a+D1)+PC.*u(y+D2)+(1-PC).*u(D2);
EUN(:,:,3) = PA.*u(x+D1)+(1-PA).*u(D1)+u(a+D2);
EUN(:,:,4) = PA.*u(x+D1)+(1-PA).*u(D1)+PC.*u(y+D2)+(1-PC).*u(D2);

DEU  = l.*(EU - repmat(max(EU,[],3),[1 1 4]));

DEUN = l.*(EUN - repmat(max(EUN,[],3),[1 1 4]));

% CPB and CPN in levels
CPB = exp(DEU)./repmat(sum(exp(DEU),3),[1 1 4]);
CPN = exp(DEUN)./repmat(sum(exp(DEUN),3),[1 1 4]);

K = repmat(broadPRi,[1 T 4]);
CP = K.*CPB + (1-K).*CPN;
CPb = CPB; CPn = CPN;

end