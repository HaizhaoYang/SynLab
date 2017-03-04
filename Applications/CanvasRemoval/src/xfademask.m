function mask = xfademask(patch_size, overlap_width, type)

switch type
    case 'black'
        mask = zeros(size(patch_size)); return;
    case 'white'
        mask = ones(size(patch_size)); return;
    case 'linear'
        horzrise = (1:(overlap_width(2)))/(overlap_width(2)+1);
        row = [horzrise, ones(1,patch_size(2)-2*overlap_width(2)), horzrise(end:-1:1)];
        % total width = patch_size(2)
        % rising    1:overlap_width(2)
        % constant  overlap_width(2)+1 : patch_size(2)-overlap_width(2)
        % falling   patch_size(2)-overlap_width(2)+1:patch_size(2)      ==end
        
        vertrise = (1:(overlap_width(1)))'/(overlap_width(1)+1);
        col = [vertrise; ones(patch_size(1)-2*overlap_width(1),1); vertrise(end:-1:1)];
        
        mask = col*row;
    otherwise
        error lolwat;
end

if patch_size(1)<2*overlap_width(1) || patch_size(2)<2*overlap_width(2)
    error lol;
end
