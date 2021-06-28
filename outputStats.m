function [M,SE,CI,STARS] = outputStats(X,sig,stars)

[R,C,S] = size(X);

for rr = 1:R
    for cc = 1:C
        v = reshape(X(rr,cc,:),[],1,1);
        M(rr,cc) = mean(v);
        SE(rr,cc) = std(v);
        CI(rr,cc,1:2) = quantile(v,[sig./2 1-sig./2]);
        s = '';
        for ss = 1:numel(stars)
            ci = quantile(v,[stars(ss)./2 1-stars(ss)./2]);
            if 0<min(ci) || 0>max(ci)
                s = [s '*'];
            end
        end
        STARS{rr,cc} = s;
    end
end


end

