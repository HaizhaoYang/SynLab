function [imgblk] = downsize(img,blocksize,overlap)

imgsize = size(img);
if (imgsize(1) <= blocksize(1) || imgsize(2) <= blocksize(2))
    imgblk = img;
    return; 
end
patch_mn = (imgsize+overlap) ./ (blocksize-overlap);
padding = round((blocksize-overlap).*(1-mod(patch_mn,1)));

thispaint = wextend(2,'sym',img,overlap+floor(padding/2),'lu');
thispaint = wextend(2,'sym',thispaint,overlap+floor(padding/2)+mod(padding,2),'rd');

paintsize = size(thispaint); 

gridplan = xfadegrid(paintsize(1:2),blocksize,overlap);

imgblk = zeros(blocksize(1),blocksize(1),length(gridplan));

for p=1:length(gridplan)
    
    thisplan = squeeze(gridplan(p,:,:));
    min_r = thisplan(1,1);      min_c = thisplan(1,2);
    max_r = thisplan(2,1);      max_c = thisplan(2,2);
    imgblk(:,:,p) = thispaint(min_r:max_r, min_c:max_c, :);
    
end
