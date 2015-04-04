function fq_out = filterFreq(fq_in,is_min)
% Filtering in the frequency domain
%
% By Haizhao Yang

[num N] = size(fq_in);
fq_out = zeros(num,N);
if ~is_min
for cnt = 1:num
  %  figure; subplot(2,2,1);plot(fq_in(cnt,:));
    temp = fftshift(fft(ifftshift(fq_in(cnt,:))))/N;
  %  subplot(2,2,2);plot(real(temp));
    temp(min(round(N/2+mean(fq_in(cnt,:))/2),N):end) = 0;
    temp(1:max(1,round(N/2-mean(fq_in(cnt,:))/2))) = 0;
  %  subplot(2,2,4);plot(real(temp));
    fq_out(cnt,:) = real(fftshift(ifft(ifftshift(temp)))*N);
    fq_out(cnt,:) = fq_out(cnt,:)*norm(fq_in(cnt,:))/norm(fq_out(cnt,:));
  %  subplot(2,2,3);plot(fq_out(cnt,:));
  %  pause;
end
else
    for cnt = 1:num
  %  figure; subplot(2,2,1);plot(fq_in(cnt,:));
    temp = fftshift(fft(ifftshift(fq_in(cnt,:))))/N;
  %  subplot(2,2,2);plot(real(temp));
    temp(min(round(N/2+min(fq_in(cnt,:))/2),N):end) = 0;
    temp(1:max(1,round(N/2-min(fq_in(cnt,:))/2))) = 0;
  %  subplot(2,2,4);plot(real(temp));
    fq_out(cnt,:) = real(fftshift(ifft(ifftshift(temp)))*N);
    fq_out(cnt,:) = fq_out(cnt,:)*norm(fq_in(cnt,:))/norm(fq_out(cnt,:));
  %  subplot(2,2,3);plot(fq_out(cnt,:));
  %  pause;
    end
end
      