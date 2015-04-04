function [dis  dis2] = distance_linear_fit(sig1,sig2)
% Measure the distance of two signals
%
% By Haizhao Yang

N = length(sig1);
x = 0:1/N:(N-1)/N;
if 0
max1 = max(abs(sig1));
max2 = max(abs(sig2));
if max1 > max2
    y = sig1./sig2;
else 
    y = sig2./sig1;
end
end
sig1 = sig1/mean(sig1);
sig2 = sig2/mean(sig2);
y = sig1./sig2;
if 0
    [p stat] = robustfit(x,y);
    dis = abs(p(2));
    yfit = polyval([p(2) p(1)],x);
else
    p = polyfit(x,y,1);
    dis = abs(p(1));
    yfit = polyval(p,x);
end
yresid = y - yfit;
dis2 = norm(yresid);



