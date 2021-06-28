function [ll] = llQREMixture2(params)

r = exp(params(1));
l = exp(params(2));
c = normcdf(params(3));

global GB G1 G2 EmpPr NT d1 d2 ActComb uid ITREATMENT

    Ub = zeros(size(EmpPr));
     Un = Ub;
    for aa = 1:4
        Ub(:,aa) = sum(EmpPr.*(repmat(GB(aa,:),[NT 1])+repmat(d1+d2,[1 4])).^r,2);
        Un(:,aa) = sum(EmpPr.*(repmat(G1(aa,:),[NT 1])+repmat(d1,[1 4])).^r,2)+sum(EmpPr.*(repmat(G2(aa,:),[NT 1])+repmat(d2,[1 4])).^r,2);
        
    end
    Ub = l.*(Ub-repmat(max(Ub,[],2),[1 4]));
    Pb = exp(Ub)./repmat(sum(exp(Ub),2),[1 4]);
    Un = l.*(Un-repmat(max(Un,[],2),[1 4]));
    Pn = exp(Un)./repmat(sum(exp(Un),2),[1 4]);
    
    llitb = sum(ActComb.*log(Pb),2);
    llitn = sum(ActComb.*log(Pn),2);
    
    uidList = unique(uid(ITREATMENT));
    lType = zeros(numel(uidList),2);
    
    
    for ii = 1:numel(uidList)
        lType(ii,1) = exp(sum(llitb(uid==uidList(ii))));
        lType(ii,2) = exp(sum(llitn(uid==uidList(ii))));
    end
    
    
    ll = sum(log(c.*lType(:,1)+(1-c).*lType(:,2)));


end

