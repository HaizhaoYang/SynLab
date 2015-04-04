function [comp amplt_cos amplt_sin] = opt_amplt(f,amplt,pre_phase)
% Solve the least square optimization in Step 6 in the diffeomorphism
% spectral analysis method.
%
% By Haizhao Yang

a = amplt.*cos(2*pi*pre_phase);
b = amplt.*sin(2*pi*pre_phase);
A = [a*transpose(a) a*transpose(b); a*transpose(b) b*transpose(b)];
d = [a*transpose(f); b*transpose(f)];
amplt_cos = (A(2,2)*d(1)-A(1,2)*d(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1));
amplt_sin = (A(2,1)*d(1)-A(1,1)*d(2))/(A(2,1)*A(1,2)-A(1,1)*A(2,2));
comp = amplt_cos.*a + amplt_sin.*b;
