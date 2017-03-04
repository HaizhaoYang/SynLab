function [fftimg_crop_cleaned, maskcrop] = checkPeak(fftimg_crop)

[ch,cw] = size(fftimg_crop);
maskcrop = zeros(size(fftimg_crop));
fftimg_crop_cleaned = fftimg_crop;
fftimg_crop = log(abs(fftimg_crop));

cy = floor(ch/2)+1;
cx = floor(cw/2)+1;

ratio = .9;
maxval = max(max(fftimg_crop));

fftimg_crop_cleaned(cy,cx) = 0;
maskcrop(cy,cx) = 1;

for i = 1:cy-1
    if fftimg_crop(cy-i,cx) > ratio*maxval && fftimg_crop(cy+i,cx) > ratio*maxval
        maskcrop(cy-i,cx) = 1;
        maskcrop(cy+i,cx) = 1;
        fftimg_crop_cleaned(cy-i,cx) = 0;
        fftimg_crop_cleaned(cy+i,cx) = 0;
    else
        break;
    end
end

for i = 1:cx-1
    if fftimg_crop(cy,cx-i) > ratio*maxval && fftimg_crop(cy,cx+i) > ratio*maxval
        maskcrop(cy,cx-i) = 1;
        maskcrop(cy,cx+i) = 1;
        fftimg_crop_cleaned(cy,cx-i) = 0;
        fftimg_crop_cleaned(cy,cx+i) = 0;
    else
        break;
    end
end


