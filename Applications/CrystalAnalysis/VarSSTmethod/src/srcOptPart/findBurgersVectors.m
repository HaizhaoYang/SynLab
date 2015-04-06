function prescribedCurl = findBurgersVectors( curlG, latticeVectors, conCompCoarse, conCompFine )
% identify Burgers vectors in coarse defect regions
curl1 = curlG(:,:,1);
curl2 = curlG(:,:,2);
defectRegions = zeros(size(curlG(:,:,1)));
BurgersVecLattice = [];
for ii = 1:conCompCoarse.NumObjects
    curlImg = [sum(curl1(conCompCoarse.PixelIdxList{ii}));sum(curl2(conCompCoarse.PixelIdxList{ii}))];
    BurgersVecLattice = [BurgersVecLattice,round(latticeVectors\curlImg)];
    defectRegions(conCompCoarse.PixelIdxList{ii}) = ii;
end

% compute total curl in fine defect regions
subComp = cell(conCompCoarse.NumObjects,1);
curlImg = [];
for ii = 1:conCompFine.NumObjects
    curlImg = [curlImg, [sum(curl1(conCompFine.PixelIdxList{ii}));sum(curl2(conCompFine.PixelIdxList{ii}))]];
    coarseComp = defectRegions(conCompFine.PixelIdxList{ii}(1));
    subComp{coarseComp} = [subComp{coarseComp},ii];
end

% distribute curl in fine defect regions so that it matches the Burgers vectors of the coarser regions
BurgersVec = curlImg;
for ii = 1:conCompCoarse.NumObjects
    curlImgLattice = latticeVectors\curlImg(:,subComp{ii});
    A = diag([1./curlImgLattice(:).^2;0;0]);
    %A = eye(2*length(subComp{ii})+2);
    A(end-1,1:2:end) = 1;
    A(end,2:2:end) = 1;
    A(1:2:end,end-1) = 1;
    A(2:2:end,end) = 1;
    A(end-1:end,end-1:end) = 0;
    b = [1./curlImgLattice,BurgersVecLattice(:,ii)];
    %b = [curlImgLattice,BurgersVecLattice(:,ii)];
    corrCurlFine = reshape(A\b(:),2,[]);
    corrCurlFine = round(corrCurlFine(:,1:end-1));
    diff = BurgersVecLattice(:,ii)-sum(corrCurlFine,2);
    for dim = 1:2
        [~,ind] = sort((curlImgLattice(dim,:)-corrCurlFine(dim,:))./curlImgLattice(dim,:));
        if diff(dim) < 0
            corrCurlFine(dim,ind(1:-diff(dim))) = corrCurlFine(dim,ind(1:-diff(dim)))-1;
        else
            corrCurlFine(dim,ind(end-diff(dim)+1:end)) = corrCurlFine(dim,ind(end-diff(dim)+1:end))+1;
        end
    end
    BurgersVec(:,subComp{ii}) = latticeVectors*corrCurlFine;
end

% compute prescribed curl
prescribedCurl1 = zeros(size(defectRegions));
prescribedCurl2 = prescribedCurl1;
for ii=1:conCompFine.NumObjects
    prescribedCurl1(conCompFine.PixelIdxList{ii}) = curl1(conCompFine.PixelIdxList{ii}) * BurgersVec(1,ii)/curlImg(1,ii);
    prescribedCurl2(conCompFine.PixelIdxList{ii}) = curl2(conCompFine.PixelIdxList{ii}) * BurgersVec(2,ii)/curlImg(2,ii);
end
prescribedCurl = cat(3,prescribedCurl1,prescribedCurl2);
end