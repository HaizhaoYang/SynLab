function img = upsize(imgblk,blocksize,overlap,imgsize)

if size(imgblk,3)==1
    img = imgblk;
    return; 
end

patch_mn = (imgsize+overlap) ./ (blocksize-overlap);
padding = round((blocksize-overlap).*(1-mod(patch_mn,1)));
paintsize = imgsize+2*overlap+2*floor(padding/2)+mod(padding,2); 

gridplan = xfadegrid(paintsize(1:2),blocksize,overlap);

patchmask = xfademask(blocksize,overlap,'linear');
img = zeros(paintsize(1),paintsize(2));

for j=length(gridplan):-1:1
    
    thisplan = squeeze(gridplan(j,:,:));
    min_r = thisplan(1,1);    min_c = thisplan(1,2);
    max_r = thisplan(2,1);    max_c = thisplan(2,2);

    thispatch = imgblk(:,:,j);
    img(min_r:max_r, min_c:max_c) = img(min_r:max_r, min_c:max_c) + thispatch.*patchmask;
    
end

img = img((overlap(1)+floor(padding(1)/2)+1):(overlap(1)+floor(padding(1)/2)+imgsize(1)), ...
    (overlap(2)+floor(padding(2)/2)+1):(overlap(2)+floor(padding(2)/2)+imgsize(2)),:);
