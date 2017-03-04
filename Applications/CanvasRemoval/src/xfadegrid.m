function [grid, npatches] = xfadegrid(mn,patch_size,overlap_width)

unit_size = patch_size - overlap_width;

patch_mn = (mn-overlap_width)./unit_size;

npatches = patch_mn(1)*patch_mn(2);
grid = zeros(npatches, 2, 2);   % upperleft(m,n), lowerright(m,n)
cursor = 1; % this is dumb
for i=1:patch_mn(1)
    for j=1:patch_mn(2)
        grid(cursor,1,:) = [1+unit_size(1)*(i-1);         1+unit_size(2)*(j-1)];    % upper-left corner
        grid(cursor,2,:) = [unit_size(1)*(i-1)+patch_size(1);    unit_size(2)*(j-1)+patch_size(2)];        % lower-right corner   
        cursor = cursor+1;
    end
end

if cursor ~= npatches+1
    disp lol;
end