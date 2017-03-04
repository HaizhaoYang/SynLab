function [res cns] = fqPeak(txt,lowFq,numPeaks,pixNum,is_display,is_stop)
% Input:
% txt - input texture image
% lowFq - low-pass parameter, do not process data in the low frequency part
% numPeaks - number of peaks need to remove
% pixNum - the radius of pixels that around a peak
% is_display - whether display results
%
% Output:
% res - residual of the input texture image that belong to the natural image
% cnx - recovered canvas from the input texture image

if nargin < 6, is_stop = 0; end;
if nargin < 5, is_display = 0; end;
if nargin < 4, pixNum = 8; end;
if nargin < 3, numPeaks = 24; end;
if nargin < 2, lowFq = 50; end;

sz = size(txt);
res = zeros(sz);
cns = zeros(sz);
if length(sz) == 2
    L = 1;
else
    L = 3;
end
for cntimg = 1:L
    txth = fftshift(fft2(ifftshift(txt(:,:,cntimg))));
    if is_display
        %figure;imagesc((txt(:,:,cntimg)));colorbar;axis image;colormap gray;
        figure;imagesc(log(abs(txth)));colorbar;axis image;
        if is_stop
            figure;
        end
    end
    
    temp = txth(round(end/2)-lowFq:round(end/2)+lowFq,round(end/2)-lowFq:round(end/2)+lowFq);
    txth(round(end/2)-lowFq:round(end/2)+lowFq,round(end/2)-lowFq:round(end/2)+lowFq) = 0;
    for cnt = 1:numPeaks
        if is_stop
            imagesc(log(abs(txth)));axis image;
        end
        maxval = max(abs(txth(:)));
        [pos1 pos2] = find(abs(txth)==maxval);
        pos1 = pos1(1); pos2 = pos2(1);
        rowNum = max(1,pos1-pixNum):min(pos1+pixNum,sz(1));
        colNum = max(1,pos2-pixNum):min(pos2+pixNum,sz(2));
        txth(rowNum,colNum) = 0;
        if is_stop
            pause;
        end
    end
    if is_display
        figure;imagesc(log(abs(txth)));axis image;
    end
    txth(round(end/2)-lowFq:round(end/2)+lowFq,round(end/2)-lowFq:round(end/2)+lowFq) = temp;
    res(:,:,cntimg) = real(fftshift(ifft2(ifftshift(txth))));
    cns(:,:,cntimg) = txt(:,:,cntimg)-res(:,:,cntimg);
    if 0
        figure;imagesc((res(:,:,cntimg)));colorbar;axis image;colormap gray;
        pic = figure;imagesc(cns(:,:,cntimg));colorbar;axis image;colormap gray;
    end
end