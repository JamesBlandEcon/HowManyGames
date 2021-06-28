function [ll] = llQRENarrow(params)

r = exp(params(1));
l = exp(params(2));

global GB G1 G2 EmpPr NT d1 d2 ActComb ITREATMENT

    U = zeros(size(EmpPr));
    for aa = 1:4
        U(:,aa) = sum(EmpPr.*(repmat(G1(aa,:),[NT 1])+repmat(d1,[1 4])).^r,2)+sum(EmpPr.*(repmat(G2(aa,:),[NT 1])+repmat(d2,[1 4])).^r,2);
    end
    U = l.*(U-repmat(max(U,[],2),[1 4]));
    
    P = exp(U)./repmat(sum(exp(U),2),[1 4]);

    ll = sum(log(P(ActComb(ITREATMENT,:)==1)));


end

