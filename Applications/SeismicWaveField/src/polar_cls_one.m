function [tb,center_pos,group] = polar_cls_one(is_real,a,b,M,NB,N,eps,lemda,R,max_group)
%This code cluster points in 2D momentum space at each location (a,b)
%is_real = 1, if the data of mode decomposition is real.
%          0, otherwise.
%M is a matrix of size NB*NB, storing the value synchrosqueezed energy
%distribution in the momentum space of size NB*NB at (a,b).
%Use NB, N, eps, and lemda to determine angle threshold theta
%R is the radius threshold
%The code would return clusters with a maximum number of group equal to
%max_group. Other information would be discarded
if nargin < 9
    max_group = 5;
end
theta = asin(eps*lemda*N/NB);
val = max(max(M));
tb = [];
[ind1 ind2] = find(M>0);
ind_R = sqrt((ind1-NB/2).^2+(ind2-NB/2).^2);
group = 0;
center_pos = [];
while group < max_group & val(1) > 0
    group = group + 1;
    [pos1 pos2] = find(M==val(1));
    center_pos = [center_pos; [pos1(1) pos2(1)]];
    
    %find neighbor in angle
    temp_agl = vec_angle(ind1-NB/2,ind2-NB/2,pos1(1)-NB/2,pos2(1)-NB/2);
    if ~is_real
        loc = find(temp_agl<=theta);
    else
        loc = find(temp_agl<theta | temp_agl>pi-theta);
    end
    %find neighbor in radius
    loc2 = find(abs(ind_R(loc)-sqrt((pos1(1)-NB/2)^2+(pos2(1)-NB/2)^2))<R);
    loc = loc(loc2);
    length_loc = length(loc);
    tb = [tb [ind1(loc)'; ind2(loc)'; ones(1,length_loc)*(pos1(1)+(pos2(1)-1)*NB)]];
    M(ind1(loc),ind2(loc)) = 0;
    val = max(max(M));
    [ind1 ind2] = find(M>0);
end
center_pos = [[a;b]*ones(1,group); center_pos'];

