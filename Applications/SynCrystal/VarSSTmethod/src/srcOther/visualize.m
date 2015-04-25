function visualize(u,G,curlG,masses,thre,binaryIndex,IfShow,address,str,example,CC)
posBD = find(binaryIndex == 0);
pic = figure;
for ii = 1:2
    subplot('Position',[(-1+ii)/6+.05,3/4,.15,.2]);
    imagesc(G(:,:,1,ii)); axis image;colorbar;
end
for ii = 1:2
    subplot('Position',[(-1+ii)/6+.05,1/2,.15,.2]);
    imagesc(G(:,:,2,ii)); axis image;colorbar;
end

% curl G
subplot('Position',[1/3+.05,1/2,.3,.4]);
imagesc(curlG(:,:,1)); axis image;colorbar;

title('curl(G)_1');
subplot('Position',[2/3+.05,1/2,.3,.4]);
imagesc(curlG(:,:,2)); axis image;colorbar;
title('curl(G)_2');



% input image
subplot('Position',[.05,.05,.3,.4]);
imagesc(u); axis image;colorbar;


% volume distortion
subplot('Position',[1/3+.05,0.05,.3,.4]);
vol = G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3);
if IfShow
    pos = find(binaryIndex==0);
    vol(pos) = 1;
end
imagesc(vol-1); axis image;  h=colorbar; caxis([-thre thre]);
colorbar;
title('Distortion Volume');


% norm curl G
subplot('Position',[2/3+.05,.05,.3,.4]);
normCurl = sqrt(sum(curlG.^2,3));
if IfShow
    normCurl(pos) = 0;
end
imagesc(normCurl); axis image; 
title('norm(curl(G))');
colorbar;


drawnow;

if 1
    if nargin > 7
        [m n] = size(G(:,:,1,1));
        bwMask = zeros(m,n);
        bwMask(posBD) = 1;
        agl = zeros(m,n);
        diff = zeros(m,n);
        for cnt = 1:m
            for cnt2 = 1:n
                temp = [G(cnt,cnt2,1,1) G(cnt,cnt2,1,2);G(cnt,cnt2,2,1) G(cnt,cnt2,2,2)];
                [U H] = polarDec(temp);
                agl(cnt,cnt2) = mod(real(acos(U(1,1))),pi/3)*180/pi;
                [ev ew] = eig(H);
                diff(cnt,cnt2) = abs(ew(1,1)-ew(2,2));
            end
        end
        agl(posBD) = 0;
        pic = figure;imagesc(agl); axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);
        colorbar;colormap(pmkmp());
        alphamask(bwMask,[0.3 0.3 0.3],1);
        title('Grain Orientation');
        
        diff(posBD) = 0;
        pic = figure;imagesc(diff); axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);colorbar;
        alphamask(bwMask,[0.3 0.3 0.3],1);
        name = sprintf('%sDiff%sEx%d.fig',address,str,example);
        if example == 2 | example == 6
            caxis([0 0.075]);%0.1
        else if example == 13
                caxis([0 0.3]);
            else if example == 1 
                    caxis([0 0.05]);
                end
            end
        end
        title('Difference in principle stretches');
        
        pic = figure;
        % volume distortion
        vol = G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3);
        if IfShow
            pos = find(binaryIndex==0);
            vol(pos) = 1;
        end
        imagesc(vol-1); axis image;  h=colorbar; caxis([-thre thre]);colorbar;
        title('Distortion Volume');
        set(gca,'xtick',[]);set(gca,'ytick',[]);
        name = sprintf('%s%sVolEx%d.fig',address,str,example);
        alphamask(bwMask,[0.3 0.3 0.3],1);
        
        
        % norm curl G
        pic = figure;
        normCurl = sqrt(sum(curlG.^2,3));
        if IfShow
            normCurl(pos) = 0;
        end
        
        if str(1) =='O'
            for i=1:CC.NumObjects
                avgCurl = sum(normCurl(CC.PixelIdxList{i}))/length(CC.PixelIdxList{i});
                normCurl(CC.PixelIdxList{i}) = avgCurl;
            end
        else
            avgCurl = sum(normCurl(posBD))/length(posBD);
            pos = find(normCurl>=avgCurl);
            normCurl(pos) = avgCurl;
        end
        imagesc(normCurl); axis image;
        title('norm(curl(G))');
        colorbar;
        set(gca,'xtick',[]);set(gca,'ytick',[]);
        colormap (1-gray);
        if example == 1
            caxis([0 0.3]);
        else if example == 2
                caxis([0 0.3]);
            else if example == 6
                    caxis([0 0.3]);
                else if example == 13
                        caxis([0 0.3]);
                    end
                end
            end
        end
        alphamask(zeros(m,n),[0.3 0.3 0.3],0);
        name = sprintf('%s%sCurlEx%d.fig',address,str,example);
        
    end
end
