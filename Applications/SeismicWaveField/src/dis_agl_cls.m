function [cls_result, num_group] = dis_agl_cls(is_real,tb,R,theta,t1,NB)
%This code cluster points in 4D phase space
%is_real = 1, if the data of mode decomposition is real.
%          0, otherwise.
%tb stores position information of points to be clustered
%R is the distance separation parameter in space, providing threshold
%distance.
%theta is the angle separation parameter in phase, providing threshold
%angle.
%t1 is the size of grid in space
%NB is the size of grid in phase
%
% By Haizhao Yang

Num_pt = size(tb,2);
group = zeros(1,Num_pt);
num_group = 0;
for cnt_pt = 1:Num_pt-1 %be careful of last pt
    loc = (cnt_pt+1):Num_pt;
    dis_vec = tb(1,loc) - tb(1,cnt_pt); bad = find(dis_vec>t1/2); dis_vec(bad) = t1-dis_vec(bad);  bad = find(dis_vec<-t1/2); dis_vec(bad) = t1+dis_vec(bad);
    loc1 = find( abs( dis_vec )<=R );
    dis_vec = tb(2,loc(loc1)) - tb(2,cnt_pt);  bad = find(dis_vec>t1/2); dis_vec(bad) = t1-dis_vec(bad);  bad = find(dis_vec<-t1/2); dis_vec(bad) = t1+dis_vec(bad);
    loc2 = find( abs( dis_vec )<=R );
    
    agl_vec = vec_angle( tb(3,loc(loc1(loc2)))-NB/2, tb(4,loc(loc1(loc2)))-NB/2, tb(3,cnt_pt)-NB/2, tb(4,cnt_pt)-NB/2);
    if ~is_real
        loc3 = find( agl_vec<=theta );%no abs
    else
        loc3 = find( agl_vec<=theta | agl_vec>=pi-theta );
    end
    vec = [group(cnt_pt) group( loc(loc1(loc2(loc3))) )];
    tmp = vec(find(vec>0));
    if length(tmp)==0
        num_group = num_group + 1;
        group( loc(loc1(loc2(loc3)) )) = num_group;
        group(cnt_pt) = num_group;
    else
        max_num = max(tmp);
        min_num = min(tmp);
        if min_num~=max_num
            for cnt_check = min_num+1:max_num
                if find( vec == cnt_check)
                    loc_check = find(group == cnt_check);
                    group(loc_check) = min_num;
                end
            end
        end
        group( loc(loc1(loc2(loc3)) )) = min_num;
        group(cnt_pt) = min_num;
    end
end
if group(end) == 0
    num_group = num_group + 1;
    group(end) = num_group;
end
gap=0;
new_group=zeros(1,Num_pt);
for i=1:max(group)
    loc=find(group==i);
    if length(loc)>0
        new_group(loc)=(i-gap)*ones(1,length(loc));
    else
        gap=gap+1;
    end
end 
group = new_group;
cls_result = [tb;group];
    