function [wl,wr] = gdct_cos(x)

% gdct_cos.m - Creates the two halves of a compactly supported window
%
% Inputs
%   x       vector or matrix of abscissae, the relevant ones from 0 to 1
%
% Outputs
%   wl,wr   vector or matrix containing samples of the left, resp. right
%           half of the window
%
% Used at least in gdct_fwd.m, gdct_ext_1.m, gdct_ext_2.m and gfdct_inv.m
%
% By Haizhao Yang

wr = zeros(size(x));
wl = zeros(size(x));
x(abs(x) < 2^-52) = 0;
wr((x > 0) & (x < 1)) = (1 + cos(pi*x((x > 0) & (x < 1))))/2  ;
wr(x <= 0) = 1;
wl((x > 0) & (x < 1)) = (1 + cos(pi*(1-x((x > 0) & (x < 1)))))/2;
wl(x >= 1) = 1;
normalization = sqrt(wl.^2 + wr.^2);
wr = wr ./ normalization;
wl = wl ./ normalization;

