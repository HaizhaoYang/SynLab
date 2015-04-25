function [CC,colComponents] = periodicConnComp( binaryImg )
CC = bwconncomp(1-binaryImg);
colComponents = zeros(size(binaryImg));
for i=1:CC.NumObjects
    colComponents(CC.PixelIdxList{i}) = i;
end
CC = ConnCompPeriodic(colComponents, CC);
for i=1:CC.NumObjects
    colComponents(CC.PixelIdxList{i}) = i;
end
end