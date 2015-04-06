function [L,G0] = findPointGroupJumpsBySweeping( pointGroup, G1, defectConComp )
% sweep all pixels via fast marching method and reorient the matrix at each pixel to fit its neighbors; first inside defect regions, then outside
G0 = G1;
queueToDo = [1,0,0,1]; % elmt idx, dist to nghbr, nghbr, ptgrpel
fixed = 0*G1(:,:,1,1);
todoOrFixed = 0*G1(:,:,1,1);
while ~isempty(queueToDo)
    %count = sum(fixed(:))
    idx = queueToDo(end,1);
    row = mod(idx-1,size(G0,1))+1;
    col = (idx-row)/size(G0,1)+1;
    ptGrp = pointGroup(:,:,queueToDo(end,4));
    queueToDo = queueToDo(1:end-1,:);
    if ~fixed(idx)
        G0(row,col,:,:) = ptGrp * squeeze(G1(row,col,:,:));
        idcs = [mod(row,size(G0,1))+1,col;...
            row,mod(col,size(G0,2))+1;...
            mod(row-2,size(G0,1))+1,col;...
            row,mod(col-2,size(G0,2))+1];
        for ii = 1:4
            if ~todoOrFixed(idcs(ii,1),idcs(ii,2))
                [ptGrp,dst] = findPointGroupConnection(squeeze(G0(row,col,:,:)),squeeze(G1(idcs(ii,1),idcs(ii,2),:,:)),pointGroup);
                if defectConComp(idcs(ii,1),idcs(ii,2))
                    dst = dst-1e3;
                elseif defectConComp(row,col)
                    dst = dst+1e3;
                end
                queueToDo = [queueToDo; [idcs(ii,1)+(idcs(ii,2)-1)*size(G0,1) dst idx ptGrp] ];
                todoOrFixed(idcs(ii,1),idcs(ii,2)) = 1;
            end
        end
        fixed(idx) = 1;
        if false %mod(count,289) == 0
            subplot(1,3,1); imagesc(fixed+2*logical(defectConComp)); axis image;
            keyboard
        end
    end
    [~,ind] = sort(queueToDo(:,2),'descend');
    queueToDo = queueToDo(ind,:);
end

% introduce jump wherever the newly aligned G has a jump (but ignore too large jumps, since those are sitting on defects and are not reliable)
L = [];
for ii = 1:size(G0,1)
    for jj = 1:size(G0,2)
        jp = mod(jj,size(G0,2))+1;
        ip = mod(ii-2,size(G0,1))+1;
        [ptGrpIdx,dst] = findPointGroupConnection( squeeze(G0(ii,jj,:,:)), squeeze(G0(ii,jp,:,:)), pointGroup );
        if ( ptGrpIdx > 1 ) && ( dst < .2 ) && ~(defectConComp(ii,jj)&&defectConComp(ii,jp))
            L = [L;[ii jj 1 ptGrpIdx]];
        end
        [ptGrpIdx,dst] = findPointGroupConnection( squeeze(G0(ii,jj,:,:)), squeeze(G0(ip,jj,:,:)), pointGroup );
        if ( ptGrpIdx > 1 ) && ( dst < .2 ) && ~(defectConComp(ii,jj)&&defectConComp(ip,jj))
            L = [L;[ii jj 0 ptGrpIdx]];
        end
    end
end
end