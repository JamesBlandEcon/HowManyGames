function [FEU,FCE,FCEabs,safeCE] = FEU_RD(rRDtrans,kRD,gam,l,ACTCOMB,ACTCOMBPm,D1,D2)

rRD = rRDtrans; % transformed r is logged

[N,T] = size(ACTCOMB);
u = @(x) x.^(repmat(rRD,[1 T]));

a = 10; c = 56; x = 100; y = 160;

% construct beliefs matrix
g = repmat(gam,[1 1 4]);
%r = repmat(rRD,[N,T,4]);

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

lEU  = l.*(EU - repmat(max(EU,[],3),[1 1 4]));
lEUN = l.*(EUN - repmat(max(EUN,[],3),[1 1 4]));

CPB = exp(lEU)./repmat((sum(exp(lEU),3)),[1 1 4]);

CPN = exp(lEUN)./repmat((sum(exp(lEUN),3)),[1 1 4]);

[temp,IB] = max(EU,[],3);
[temp,IN] = max(EUN,[],3);

FEU = zeros(N,T);
FCE = zeros(N,T);
FCEabs = zeros(N,T);
safeCE = zeros(N,T);

for tt = 1:T
    for ii = 1:N
        FEU(ii,tt) = (EU(ii,tt,IB(ii,tt))-EU(ii,tt,IN(ii,tt)));
        
        FCE(ii,tt) = 1-(EU(ii,tt,IN(ii,tt))./EU(ii,tt,IB(ii,tt))).^(1./rRD(ii));
        FCEabs(ii,tt) = EU(ii,tt,IB(ii,tt)).^(1./rRD(ii))-EU(ii,tt,IN(ii,tt)).^(1./rRD(ii));
        safeCE(ii,tt) = 1-(EU(ii,tt,1)./EU(ii,tt,IB(ii,tt))).^(1./rRD(ii));
    end
end



end