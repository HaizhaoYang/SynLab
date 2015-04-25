function CC = ConnCompPeriodic(colComponents, CCold)

[m, n] = size(colComponents);

AdjMat = eye(CCold.NumObjects);

for i = 1:m
    if (colComponents(i, 1) > 0) && (colComponents(i, n) > 0)
        AdjMat(colComponents(i, 1), colComponents(i, n)) = 1;
        AdjMat(colComponents(i, n), colComponents(i, 1)) = 1;
    end
end

for i = 1:n
    if (colComponents(1, i) > 0) && (colComponents(m, i) > 0)
        AdjMat(colComponents(1, i), colComponents(m, i)) = 1;
        AdjMat(colComponents(m, i), colComponents(1, i)) = 1;
    end
end

[group, Num_group] = DistanceClustering(AdjMat, 0.5, CCold.NumObjects, ...
                                        'large');

CC.Connectivity = CCold.Connectivity; 
CC.ImageSize = CCold.ImageSize; 
CC.NumObjects = Num_group; 
CC.PixelIdxList = cell(1, Num_group);
for i = 1:Num_group
    CC.PixelIdxList{i} = []; 
    for j = 1:CCold.NumObjects
        if (group(j) == i)
            CC.PixelIdxList{i} = [CC.PixelIdxList{i}; CCold.PixelIdxList{j}];
        end
    end
end