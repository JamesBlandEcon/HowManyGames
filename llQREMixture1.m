function [ll] = llQREMixture1(params)

r = exp(params(1));
l = exp(params(2));
c = normcdf(params(3));

global GB G1 G2 EmpPr NT d1 d2 ActComb ITREATMENT

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
    
    P = c.*Pb+(1-c).*Pn;

    ll = sum(log(P(ActComb(ITREATMENT,:)==1)));

end

