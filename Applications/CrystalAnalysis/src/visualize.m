function visualize(u,G,curlG,masses,thre,binaryIndex,IfShow,address,str,example)
% function visualize(u,G,curlG)
pic = figure;
for ii = 1:2
    subplot('Position',[(-1+ii)/6+.05,3/4,.15,.2]);
    imagesc(G(:,:,1,ii)); axis image;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
end
for ii = 1:2
    subplot('Position',[(-1+ii)/6+.05,1/2,.15,.2]);
    imagesc(G(:,:,2,ii)); axis image;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
end

% curl G
subplot('Position',[1/3+.05,1/2,.3,.4]);
imagesc(curlG(:,:,1)); axis image;
%freezeColors;
colorbar;%cbfreeze(colorbar);

title('curl(G)_1');
subplot('Position',[2/3+.05,1/2,.3,.4]);
imagesc(curlG(:,:,2)); axis image;
%freezeColors;
colorbar;%cbfreeze(colorbar);
title('curl(G)_2');



% input image
subplot('Position',[.05,.05,.3,.4]);
imagesc(u); axis image;
%freezeColors;
colorbar;%cbfreeze(colorbar);


% volume distortion
subplot('Position',[1/3+.05,0.05,.3,.4]);
vol = G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3);
if IfShow
    pos = find(binaryIndex==0);
    vol(pos) = 1;
end
imagesc(vol-1); axis image;  h=colorbar; caxis([-thre thre]);
%freezeColors;
colorbar;%cbfreeze(colorbar);
%cbfreeze(h);
title('Distortion Volume');

%freezeColors;


% norm curl G
subplot('Position',[2/3+.05,.05,.3,.4]);
normCurl = sqrt(sum(curlG.^2,3));
if IfShow
    normCurl(pos) = 0;
end
imagesc(normCurl); axis image; %colormap(1-gray); caxis([0 .25]);
title('norm(curl(G))');
%freezeColors;
colorbar;%cbfreeze(colorbar);

%val = max(temp(:));


drawnow;


if nargin > 7
    pic = figure;
    for ii = 1:2
        subplot(2,2,ii);
        imagesc(G(:,:,1,ii)); axis image;
        %freezeColors;
        colorbar;%cbfreeze(colorbar);
    end
    for ii = 1:2
        subplot(2,2,2+ii);
        imagesc(G(:,:,2,ii)); axis image;
        %freezeColors;
        colorbar;%cbfreeze(colorbar);
    end
    name = sprintf('%sG%sEx%d.fig',address,str,example);
    saveas(pic,name);
    
    pic = figure;
    % curl G
    subplot(2,2,1);
    imagesc(curlG(:,:,1)); axis image;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
    
    title('curl(G)_1');
    subplot(2,2,2);
    imagesc(curlG(:,:,2)); axis image;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
    title('curl(G)_2');
    
    % volume distortion
    subplot(2,2,3);
    vol = G(:,:,1).*G(:,:,4)-G(:,:,2).*G(:,:,3);
    if IfShow
        pos = find(binaryIndex==0);
        vol(pos) = 1;
    end
    imagesc(vol-1); axis image;  h=colorbar; caxis([-thre thre]);
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
    %cbfreeze(h);
    title('Distortion Volume');
    
    %freezeColors;
    
    
    % norm curl G
    subplot(2,2,4);
    normCurl = sqrt(sum(curlG.^2,3));
    if IfShow
        normCurl(pos) = 0;
    end
    imagesc(normCurl); axis image; %colormap(1-gray); caxis([0 .25]);
    title('norm(curl(G))');
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
    
    name = sprintf('%s%sEx%d.fig',address,str,example);
    saveas(pic,name);
    
    % input image
    pic = figure;
    imagesc(u); axis image;
    %freezeColors;
    colorbar;%cbfreeze(colorbar);
    name = sprintf('%sImgEx%d.fig',address,example);
    saveas(pic,name);
    
end
end
