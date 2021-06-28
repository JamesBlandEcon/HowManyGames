function [ll,llwide,ll1,ll2B,A1,A2N,A2B] = LTlogit(rtrans,k,l,i1,i2,c1,c2,points)

r = exp(rtrans); % transformed r is logged

X0 = linspace(0,100,points);

NC = numel(c1);
N = numel(r);

K = repmat(k,[1 NC]);

R = repmat(r,[1 NC]);
C1 = repmat(c1,[N 1]);
C2 = repmat(c2,[N 1]);


UX1  = l.*(0.5.*(C1.*i1+100-i1).^R + 0.5.*(100-i1).^R);
UX2N = l.*(0.5.*(C2.*i2+100-i2).^R + 0.5.*(100-i2).^R);
UX2B = l.*0.25.*((2*C2.*i2+200-2.*i2).^R+2.*(C2.*i2+200-2.*i2).^R+(200-2.*i2).^R);

A1  = zeros(N,NC);
A2N = zeros(N,NC);
A2B = zeros(N,NC);

for xx = 1:(numel(X0)-1)
    xlb = X0(xx); xub = X0(xx+1);
    dx = xub-xlb;
    A1  = A1 + 0.5.*dx.*(exp(l.*(0.5.*(C1.*xub+100-xub).^R + 0.5.*(100-xub).^R))+exp(l.*(0.5.*(C1.*xlb+100-xlb).^R + 0.5.*(100-xlb).^R)));
    A2N = A2N + 0.5.*dx.*(exp(l.*(0.5.*(C2.*xub+100-xub).^R + 0.5.*(100-xub).^R))+exp(l.*(0.5.*(C2.*xlb+100-xlb).^R + 0.5.*(100-xlb).^R)));
    A2B = A2B + 0.5.*dx.*(...
                +exp(l.*0.25.*((2*C2.*xub+200-2*xub).^R+2.*(C2.*xub+200-2*xub).^R+(200-2*xub).^R))...
                +exp(l.*0.25.*((2*C2.*xlb+200-2*xlb).^R+2.*(C2.*xlb+200-2*xlb).^R+(200-2*xlb).^R))...
        );
end

% log likelihoods

ll1  = UX1-log(A1);
ll2N = UX2N-log(A2N);
ll2B = UX2B-log(A2B);

llwide = (ll1+K.*ll2B+(1-K).*(ll2N));
ll = sum(llwide,2);

end


% function [ll,llwide,ll1,ll2B,Const1,Const2N,Const2B] = LTlogit(r,k,l,i1,i2,c1,c2,order)
% 
% X0 = 10:10:90;
% X0lb = X0-10;
% X0ub = X0+10;
% 
% NC = numel(c1);
% N = numel(r);
% 
% K = repmat(k,[1 NC]);
% 
% R = repmat(r,[1 NC]);
% C1 = repmat(c1,[N 1]);
% C2 = repmat(c2,[N 1]);
% 
% 
% UX1  = l.*(0.5.*(C1.*i1+100-i1).^R + 0.5.*(100-i1).^R);
% UX2N = l.*(0.5.*(C2.*i2+100-i2).^R + 0.5.*(100-i2).^R);
% UX2B = l.*0.25.*((2*C2.*i2+200-i2).^R+2.*(C2.*i2+200-i2).^R+(200-i2).^R);
% 
% % set up taylor polynomial
% % common components
% B(:,:,1) =1;
% for oo = 2:order;
%         B(:,:,oo) = B(:,:,oo-1).*l./(oo-1).*(R+2-oo);
% end
% DU1 = zeros(N,NC,order);
% DU2N = zeros(N,NC,order);
% DU2B = zeros(N,NC,order);
% Poly0 = zeros(N,NC,order);
% 
% 
% 
% llwide = zeros(N,2*NC);
% for xx = 1:numel(X0)
%     x0 = X0(xx);
%     UX1o  = l.*(0.5.*(C1.*x0+100-x0).^R + 0.5.*(100-x0).^R);
%     UX2No = l.*(0.5.*(C2.*x0+100-x0).^R + 0.5.*(100-x0).^R);
%     UX2Bo = l.*0.25.*((2*C2.*x0+200-x0).^R+2.*(C2.*x0+200-x0).^R+(200-x0).^R);
% for oo = 1:order;
%     
% 
%     
%     
%     % derivatives of utility (without the r multiplier)
%     DU1(:,:,oo) =  UX1o.*0.5.*((C1-1).^(oo-1).*(C1.*x0+100-x0).^(R+1-oo)+(-1).^oo.*(100-x0).^(R+1-oo));
%     DU2N(:,:,oo) =  UX2No.*0.5.*((C2-1).^(oo-1).*(C2.*x0+100-x0).^(R+1-oo)+(-1).^oo.*(100-x0).^(R+1-oo));
%     DU2B(:,:,oo) =  UX2Bo.*0.25.*(...
%                     (2.*C2-2).^(oo-1).*(2*C2.*x0+200-2.*x0).^(R+1-oo)...
%                     +2.*(C2-2).^(oo-1).*(C2.*x0+200-2.*x0).^(R+1-oo)...
%                     +(-2).^(oo-1).*(200-2.*x0).^(R+1-oo));
%     % polynomial component
%     Poly0(:,:,oo) = ((X0ub(xx)-x0).^oo-(X0lb(xx)-x0).^oo)./oo;
%     
% end
% 
% end
% 
% % log likelihoods
% Const1 = sum(B.*DU1.*Poly0,3)
% Const2N = sum(B.*DU2N.*Poly0,3)
% Const2B = sum(B.*DU2B.*Poly0,3)
% 
% ll1  = UX1-log(Const1);
% ll2N = UX2N-log(Const2N);
% ll2B = UX2B-log(Const2B);
% 
% llwide = llwide + (ll1+K.*ll2B+(1-K).*(ll2N));
% ll = sum(llwide,2);
% 
% end

% function [ll,llwide] = LTlogit(r,k,l,i1,i2,c1,c2,order)
% 
% X0 = 10:10:90;
% X0lb = X0-10;
% X0ub = X0+10;
% 
% NC = numel(c1);
% N = numel(r);
% 
% K = repmat(k,[1 NC]);
% 
% R = repmat(r,[1 NC]);
% C1 = repmat(c1,[N 1]);
% C2 = repmat(c2,[N 1]);
% 
% 
% UX1 = exp(l.*(0.5.*(C1.*i1+100-i1).^(1-R)./(1-R) + 0.5.*(100-i1).^(1-R)./(1-R)));
% 
% % set up taylor polynomial
% % common components
% A(:,:,1) = 1./(1-R);
% for oo = 2:order;
%         A(:,:,oo) = A(:,:,oo-1).*l./(oo-1).*(3-R-oo);
% end
% DU1 = zeros(N,NC,order);
% DU2N = zeros(N,NC,order);
% DU2B = zeros(N,NC,order);
% Poly0 = zeros(N,NC,order);
% PolyX1 = zeros(N,NC,order);
% PolyX2 = zeros(N,NC,order);
% 
% 
% llwide = zeros(N,2*NC);
% for xx = 1:numel(X0)
%     x0 = X0(xx);
% for oo = 1:order;
%     
%     
%     % derivatives of utility (without the r multiplier)
%     DU1(:,:,oo) =  0.5.*((C1-1).^(oo-1).*(C1.*x0+100-x0).^(2-oo-R)+(-1).^oo.*(100-x0).^(2-oo-R));
%     DU2N(:,:,oo) =  0.5.*((C2-1).^(oo-1).*(C2.*x0+100-x0).^(2-oo-R)+(-1).^oo.*(100-x0).^(2-oo-R));
%     DU2B(:,:,oo) =  0.25.*(...
%                     (2.*C2-2).^(oo-1).*(2*C2.*x0+200-2.*x0).^(2-oo-R)...
%                     +2.*(C2-2).^(oo-1).*(C2.*x0+200-2.*x0).^(2-oo-R)...
%                     +(-2).^oo.*(200-2.*x0).^(2-oo-R));
%     % polynomial component
%     Poly0(:,:,oo) = ((100-x0).^oo-(0-x0).^oo)./oo;
%     PolyX1(:,:,oo) = (i1-x0).^(oo-1);
%     PolyX2(:,:,oo) = (i2-x0).^(oo-1);
%     
% end
% % log likelihoods
% ll1  = log(sum(A.*DU1.*PolyX1,3))-log(sum(A.*DU1.*Poly0,3));
% ll2N = log(sum(A.*DU2N.*PolyX2,3))-log(sum(A.*DU2N.*Poly0,3));
% ll2B = log(sum(A.*DU2B.*PolyX2,3))-log(sum(A.*DU2B.*Poly0,3));
%     llwide = llwide + (ll1+K.*ll2B+(1-K).*(ll2N))./numel(X0);
% end
% ll = sum(llwide,2);
% 
% end


%function [ll] = LTlogit(r,k,l,i1interp,i2interp,xgrid,c1,c2)
%xgs = numel(xgrid);
%StepSize = 100./(xgs-1);
%R = repmat(r,[1 4 xgs]);
%C1 = repmat(c1,[size(i1interp,1) 1 xgs]);
%C2 = repmat(c2,[size(i1interp,1) 1 xgs]);
%
%
%
%
%xg = repmat(reshape(xgrid,1,1,[]),[size(i1interp(:,:,1)) 1]);
%
%Px1 = l.*((C1.*xg + 100 - xg+0.1).^(1-R)./(1-R) + (100 - xg+0.1).^(1-R)./(1-R));
%Px1 = exp(Px1 - repmat(max(Px1,[],3),[1 1 xgs]));
%Px1 = Px1./repmat(sum(Px1,3),[1 1 xgs]);
%A1 = (repmat(sum(Px1,3) -0.5.*Px1(:,:,1)-0.5.*Px1(:,:,end),[1 1 xgs])).*StepSize;
%Px1 = Px1./A1;
%
%lb1 = sum(Px1.*(xg==repmat(i1interp(:,:,1),[1 1 xgs])),3);
%ub1 = sum(Px1.*(xg==repmat(i1interp(:,:,2),[1 1 xgs])),3);
%c1  = sum(Px1.*(xg==repmat(i1interp(:,:,3),[1 1 xgs])),3);
%ll1 = sum(log(lb1 + (ub1-lb1).*c1),2);
%
%ll = ll1;
%
%end

