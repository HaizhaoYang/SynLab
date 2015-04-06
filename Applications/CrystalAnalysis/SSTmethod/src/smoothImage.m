function v = smoothImage( u, pixelSize, stdDeviation )
myfilter = fspecial('gaussian',pixelSize * [1 1],stdDeviation);
v = imfilter( u, myfilter, 'circular' );
end