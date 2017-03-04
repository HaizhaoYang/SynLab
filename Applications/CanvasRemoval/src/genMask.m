close all;
temp = ifftshift(fft2(fftshift(PeakSig1)));
figure;imagesc(abs(temp));
temp(235:280,235:280) = 0;
PeakSig2 = abs(fftshift(ifft2(ifftshift(temp))));
figure;imagesc(PeakSig1);
figure;imagesc(PeakSig1-PeakSig2);
