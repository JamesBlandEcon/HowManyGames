function [X] = LTinterp(r,xgrid,rgrid)

X = zeros(size(r,1),4);

for cc = 1:4
    X(:,cc) = interp1(rgrid(:,1),xgrid,r,'linear','extrap');
end

end

