function [ptGrpIdx,minVal] = findPointGroupConnection( G1, G2, pointGroup )
% returns index of point group element R (in pointGroup) such that |G1 - R G2| is smallest
minVal = 1e10;
for k = 1:size(pointGroup,3)
    val = norm(G1-pointGroup(:,:,k)*G2,2);
    if val < minVal
        minVal = val;
        ptGrpIdx = k;
    end
end
end