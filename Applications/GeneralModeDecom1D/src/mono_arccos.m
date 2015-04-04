function y = mono_arccos(x)
% Compute monotonous acos(x)
%
% By Haizhao Yang

N = length(x);
y = zeros(size(x));
y(1) = acos(x(1));
k = 1;
kk = 0;
for cnt = 2:N
    temp_val = kk*2*pi + k*acos(x(cnt));
    if temp_val < y(cnt-1)
        kk = kk + mod((k+1)/2,2);
        k = -k;
        temp_val = kk*2*pi + k*acos(x(cnt));
    end
    y(cnt) = temp_val;
end