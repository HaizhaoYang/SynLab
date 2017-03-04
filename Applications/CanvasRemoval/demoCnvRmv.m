% This is a three-step method for canvas removal
%
% By 

close all;
clear;

%% Choose examples and parameters
% lowFq: low frequency components that will be ignored
% r1,r2: control the width of ridges that will be removed
% prom: prominence of peaks when looking for periodic components, this
%   should be made dependent of the peaky signal
% maxind: maximum number of ridges in horizontal and vertical direction
% th: if set to 1, multiscale thresholding will take place in frequency domain
% p: controls the thresholding

example = 1;
switch example
    case 1
        [dat map] = imread('bird_Original.tif');
        dat = double(dat);
        sigma = 6;
        blocksize = [1024,1024];
        overlap = [64,64];
        
        lowFq = 15;
        r1 = 2; r2 = 2;
        prom = 20;
        maxind = 4;
        th = 1; p = .9;
        
        method = 4;
        if method == 1
            prom = .4;
        else
            prom = 2E-2;
        end

end

sz = size(dat);

if length(sz)>2
    n = sz(3);
else
    n = 1;
end

%% Step 1: Cartoon-texture decomposition
% Input:  dat - raw X-ray image data
% Output: img - cartoon part of the natural image
%         txt - texture part of the X-ray data
% src = img + txt

[src,img,txt] = cartoonTexture(dat,sigma);

is_display = 0;

if is_display
    if length(sz) == 2
        pic = figure; imagesc(src);axis image; axis off;  colormap gray;
        head = sprintf('Example %d, Input image, sigma = %f',example,sigma); title(head);
        pic = figure; imagesc(img);axis image; axis off;  colormap gray;
        head = sprintf('Example %d, Cartoon part, sigma = %f',example,sigma);title(head);
        pic = figure; imagesc(txt);axis image; axis off;  colormap gray;
        head = sprintf('Example %d, Texture part, sigma = %f',example,sigma);title(head);
    else
        pic = figure; imshow(uint8(src));axis image; axis off;
        head = sprintf('Example %d, Input image, sigma = %f',example,sigma); title(head);
        pic = figure; imshow(uint8(img));axis image; axis off;
        head = sprintf('Example %d, Cartoon part, sigma = %f',example,sigma);title(head);
        pic = figure; imshow(double(rgb2gray(txt)));axis image; axis off;  colormap gray;
        head = sprintf('Example %d, Texture part, sigma = %f',example,sigma);title(head);
    end
end


%% Step 2: Thresholding in the frequency domain using FFT
% Input:  txt
% Output: cnv - canvas
%         res - residual of the natural image
% txt = cnv + res

is_display = 0;

if n > 1
    txtgray = rgb2gray(txt);
    txtgrayblk = downsize(txtgray,blocksize,overlap);
else
    txtgrayblk = downsize(txt,blocksize,overlap); 
end

for channel = 1:n
    txtblk = downsize(txt(:,:,channel),blocksize,overlap);
    for idx = 1:size(txtblk,3)
        % clean spectrum
        peakloc = locatePeaks(txtgrayblk(:,:,idx),prom,blocksize,maxind,method,is_display);
        [fftimg_cleaned,mask] = removePeaks(txtblk(:,:,idx),peakloc,lowFq,r1,r2,is_display,p,method);
        restxt = real(ifft2(ifftshift(fftimg_cleaned)));
        txtblk(:,:,idx) = restxt;
    end
    txt(:,:,channel) = upsize(txtblk,blocksize,overlap,sz(1:2));
end

%% Step 3: post-processing for X-rays
% bright areas
datsmooth = im2double(uint8(img));
datbw = im2bw(datsmooth,0.9);
rec = img+txt;
reccor = rec;
reccor(datbw == 1) = dat(datbw == 1);

% dark areas
datblack = 1-im2double(uint8(dat));
datbw = im2bw(datblack,0.9);
cc = bwconncomp(datbw);
numPixels = cellfun(@numel,cc.PixelIdxList);
% [ind] = find(numPixels < 100);
[ind] = find(numPixels < 1000);
for i=1:length(ind)
    datbw(cc.PixelIdxList{ind(i)}) = 0;
end 
reccor(datbw == 1) = dat(datbw == 1);

%% output
figure, imshow(uint8(dat)); title('Original')
figure, imshow(uint8(rec)); title('Reconstructed');
figure, imshow((dat-reccor),[]); colormap gray; title('Difference')