function [X] = LTinterpBoth(r,k,xgrid,rgrid1,rgrid2N,rgrid2B)

X1 = zeros(size(r,1),4);
X2 = zeros(size(r,1),4);


for cc = 1:4
    X1(:,cc) = interp1(rgrid1(:,cc),xgrid,r,'linear','extrap');
    X2(k==0,cc) = interp1(rgrid2N(:,cc),xgrid,r(k==0),'linear','extrap');
    X2(k==1,cc) = interp1(rgrid2B(:,cc),xgrid,r(k==1),'linear','extrap');
end

X = [X1 X2];

end

