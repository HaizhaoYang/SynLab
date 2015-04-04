function y = filt_hard_threshold(x,pct,is_real)
% Hard thresholding in the frequency domain
%
% By Haizhao Yang

energy = norm(x);
y=fftshift(fft(ifftshift(x)))/length(x);
loc = find(abs(y)<max(abs(y))*pct);
y(loc)=0;
if is_real
    y =real(fftshift(ifft(ifftshift(y)))*length(y));
else
    y =fftshift(ifft(ifftshift(y)))*length(y);
end
y = y*energy/norm(y);
    