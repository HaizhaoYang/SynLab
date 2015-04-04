function y = gauss(x,sigma)
% Generate a Gaussian function

if size(x,1)>1
    y = exp(-sum(x.^2)/(2*sigma^2));
else
    y = exp(-x.^2/(2*sigma^2));
end