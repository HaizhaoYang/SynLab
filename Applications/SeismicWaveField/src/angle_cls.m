function [tb,center_pos,group] = angle_cls(is_real,a,b,M,NB,N,eps,lemda,max_group)
%This code cluster points in 2D momentum space at each location (a,b)
%is_real = 1, if the data of mode decomposition is real.
%          0, otherwise.
%M is a matrix of size NB*NB, storing the value synchrosqueezed energy
%distribution in the momentum space of size NB*NB at (a,b).
%Use NB, N, eps, and lemda to determine angle threshold theta
%The code would return clusters with a maximum number of group equal to
%max_group. Other information would be assigned to the nearest cluster.
%
% By Haizhao Yang

if nargin < 9
    max_group = 5;
end
theta = asin(eps*lemda*N/NB);
val = max(max(M));
tb = [];
[ind1 ind2] = find(M>0);
group = 0;
center_pos = [];
while val(1)>0
    group = group + 1;
    if group > max_group
        break;
    end
    [pos1 pos2] = find(M==val(1));
    center_pos = [center_pos; [pos1(1) pos2(1)]];
    temp_agl = vec_angle(ind1-NB/2,ind2-NB/2,pos1(1)-NB/2,pos2(1)-NB/2);
    if ~is_real
        loc = find(temp_agl<=theta);
    else
        loc = find(temp_agl<theta | temp_agl>pi-theta);
    end
    length_loc = length(loc);
    tb = [tb [ind1(loc)'; ind2(loc)'; ones(1,length_loc)*(pos1(1)+(pos2(1)-1)*NB)]];
    M(ind1(loc),ind2(loc)) = 0;
    val = max(max(M));
    [ind1 ind2] = find(M>0);
end
if group > max_group
  %  fprintf('group is larger than max_group \n');
    group = max_group;
    [loc1 loc2] = find(M>0);
    for cnt = 1:length(loc1)
       temp_agl = min(vec_angle(center_pos(:,1)-NB/2,center_pos(:,2)-NB/2,loc1(cnt)-NB/2,loc2(cnt)-NB/2),vec_angle(-center_pos(:,1)+NB/2,-center_pos(:,2)+NB/2,loc1(cnt)-NB/2,loc2(cnt)-NB/2));
       [val pos] = min(temp_agl);
       tb = [tb [loc1(cnt);loc2(cnt); (center_pos(pos,1)+NB*(center_pos(pos,2)-1))]];
    end
end
center_pos = [[a;b]*ones(1,group); center_pos'];

