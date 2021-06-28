%This function generates draws from a truncated normal distribution
%using the inversion method. 
%
%The syntax is  
%             truncnorm3(mu,variance,a,b) 
%where a and b are the lower and upper (respectively) limits of truncation. 
%
%NOTE that mu, variance, a and b can be vectors of equal length
%(this will save time to avoid looping)
%When there is NO UPPER BOUND (i.e., b = infinity), SET b = 999 
%When there is NO LOWER BOUND (i.e., a = -infinity), SET a = -999
%
%Thus, truncnorm3(0,1,-999,999) would give a draw from the 
%standard normal distribution. 


function [draws] = truncnorm3(mu,variance,a,b);
stderrs = sqrt(variance);

points_a = find(a==-999);
points_b = find(b==999);

a_term = normcdf( (a-mu)./stderrs);
a_term(points_a) = 0;

b_term = normcdf( (b-mu)./stderrs);
b_term(points_b) = 1;

uniforms = rand(length(mu),1);

p = a_term + uniforms.*(b_term - a_term);

draws = mu + stderrs.*norminv(p);
