function [fftimg_cleaned,mask] = removePeaks(img,peakloc,lowFq,r1,r2,is_display,p,method)

N = size(img,1);
%% generate S_r1 and S_r2
peakv = unique(peakloc(:,1));% TODO: check
peakh = unique(peakloc(:,2));
S_r1 = zeros(size(img));
S_r2 = zeros(size(img));
% use r1 to create S_r1 in our paper
for cnt = 1:length(peakv)
    S_r1(peakv(cnt)-r1:peakv(cnt)+r1,:) = 1;
    S_r2(peakv(cnt)-r2:peakv(cnt)+r2,:) = 1;
end
for cnt = 1:length(peakh)
    S_r1(:,peakh(cnt)-r1:peakh(cnt)+r1) = 1;
    S_r2(:,peakh(cnt)-r2:peakh(cnt)+r2) = 1;
end

%% multiscale paritioning using [Nv Nh]
% horizontal; compute Nh using the algorithm in DSA
peak = peakh-N/2;
pos = find(peak>0);
peak = peak(pos);
if peak(1)< lowFq
    peak = peak(2:end);
end
is_mean = 0;
[Nh spech] = instant_freq_est(peak,lowFq,is_mean);
Nh = round(Nh);
% vertical; compute Nv using the algorithm in DSA
peak = peakv-N/2;
pos = find(peak>0);
peak = peak(pos);
if peak(1)< lowFq
    peak = peak(2:end);
end
is_mean = 0;
[Nv specv] = instant_freq_est(peak,lowFq,is_mean);
Nv = round(Nv);
% compute partition indices
idxv = [0:Nv:N/2 -Nv:-Nv:-N/2] + N/2 - round(Nv/2); % start vertical position of the block of size Nv*Nh
idxh = [0:Nh:N/2 -Nh:-Nh:-N/2] + N/2 - round(Nh/2); % start horizontal position of the block of size Nv*Nh

%% multiscale thresholding
imgh = fftshift(fft2(img));
imghSmoothed = smoothImage(imgh,r2,r2/2);
imghOrg = imgh;
imgRes = imgh;
mask = zeros(size(img));
is_showBlock = 0;
if is_showBlock
    figure;
end
for cntv = 1:length(idxv)
    for cnth = 1:length(idxh)
        stv = max(idxv(cntv),1); edv = min(stv + Nv,N);
        sth = max(idxh(cnth),1); edh = min(sth + Nh,N);
        % B_nm is located in [stv:edv,sth:edh]
        %[stv edv sth edh]
        S_r1Block = S_r1(stv:edv,sth:edh);
        S_r2Block = S_r2(stv:edv,sth:edh);
        imghBlock = imgh(stv:edv,sth:edh);
        maskBlock = zeros(size(imghBlock));
        pos = find(S_r2Block==1); % B_nm \cap S_r2
        if numel(pos) > 0
            if is_showBlock
                matBlk = S_r1;
                matBlk(stv:edv,sth:edh) = 2;
                imagesc(matBlk);
                pause;
            end
            pos = find(S_r2Block==0); % B_nm \setminus S_r2
            switch method
                case 1
                    val = imghSmoothed(stv:edv,sth:edh);
                    threshold = max(abs(val(pos))); % Need to optimize this step!!
                case 2
                    val = imghSmoothed(stv:edv,sth:edh);
                    threshold = max(abs(val(pos)))/2; % Need to optimize this step!!
                case 3
                    val = imgh(stv:edv,sth:edh);
                    threshold = median(abs(val(pos))); % Need to optimize this step!!
                case 4
                    val = imgh(stv:edv,sth:edh);
                    threshold = p*numel(pos);
            end
            pos = find(S_r1Block==1); % B_nm \cap S_r1
            posRemove = find(abs(imghBlock(pos))>threshold);
            imghBlock(pos(posRemove)) = 0;
            maskBlock(pos(posRemove)) = 1;
            imgRes(stv:edv,sth:edh) = imghBlock;
            mask(stv:edv,sth:edh) = maskBlock;
        end
    end
end
imgRes(N/2-lowFq:N/2+lowFq,N/2-lowFq:N/2+lowFq) = imghOrg(N/2-lowFq:N/2+lowFq,N/2-lowFq:N/2+lowFq);
%% output results
fftimg_cleaned = imgRes;

%% visualize results
if (is_display)
    figure;subplot(1,2,1);
    imagesc(S_r1);axis image;colorbar;title('S_r');
    temp = abs(imgh);
    pos = find(S_r1==1);
    temp(pos) = 0;
    temp = temp(end/4:3*end/4,end/4:3*end/4);
    val = max(temp(:));
    subplot(1,2,2);imagesc(temp);axis image;colorbar;title('spectrum outside S_r');colormap (1-gray);caxis([0 val/2]);
    figure;imagesc(temp);axis image;title('spectrum outside S_r');colormap (1-gray);caxis([0 val/2]);
    figure;
    temp = abs(imgh);
    pos = find(mask==1);
    temp(pos) = 0;
    temp = temp(end/4:3*end/4,end/4:3*end/4);
    val = max(temp(:));
    subplot(1,2,1);imagesc(mask);axis image;colorbar;title('mask');
    subplot(1,2,2);imagesc(temp);axis image;colorbar;title('spectrum outside mask');colormap (1-gray);caxis([0 val/2]);
    temp = abs(imgh(end/4:3*end/4,end/4:3*end/4)-imgRes(end/4:3*end/4,end/4:3*end/4));
    val = max(temp(:));
    figure;imagesc(temp);axis image;colorbar;title('estimated spectrum of the canvas');colormap (1-gray);caxis([0 val/2]);
    temp = abs(imgh(end/4:3*end/4,end/4:3*end/4));
    val = max(temp(:));
    figure;imagesc(temp);axis image;colorbar;title('original spectrum');colormap (1-gray);caxis([0 val/2]);
    figure;imagesc(mask(end/4:3*end/4,end/4:3*end/4));axis image;colorbar;title('mask');
end
