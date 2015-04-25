function saveImg( img, range, cmap, mask, cmask, filename )

  % rescale image to range [0,1]
  img = ( double(img) - range(1) ) / ( range(2) - range(1) );
  img = max(img,0);
  img = min(img,1);

  % translate img into RGB
  len = length(cmap);
  img = uint8(img*(len-1)+1);
  RGB{1} = cmap(img,1);
  RGB{2} = cmap(img,2);
  RGB{3} = cmap(img,3);

  % add mask
  for i = 1:3
    RGB{i}(mask) = cmask(i);
  end
  
  % save image
  img = reshape([RGB{1},RGB{2},RGB{3}],[size(img),3]);
  imwrite(img,filename);

end