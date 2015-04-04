function amplt = amplt_est(is_real,y)
% Estimate instantaneous frequency by computing analytic signals
%
% By Haizhao Yang

if is_real
    yh = fftshift(fft(ifftshift(y)));
    yh(1:floor(end/2)) = 0;
    yh(floor(end/2)+1:end) = 2*yh(floor(end/2)+1:end);
    y = fftshift(ifft2(ifftshift(yh)));
    amplt = abs(y);
else
    amplt = abs(y);
end