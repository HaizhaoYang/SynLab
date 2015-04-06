function [group,Num_group] = DistanceClustering(W,control,goal_num,option)
Num_pt = size(W,2);
if Num_pt == 1
    Num_group = 1;
    group = 1;
else
    W = W+eye(Num_pt)*control;
    group = zeros(Num_pt,1);
    Num_group = 1;
    for cnt = 1:Num_pt
        if option == 'small'
           loc = find(W(cnt,:)<=control);
        else
            loc = find(W(cnt,:)>=control);
        end
        [minval,minpos] = min(group(loc));
        [maxval,maxpos] = max(group(loc));
        if maxval == 0
           group(loc) = Num_group*ones(1,length(loc));
           Num_group = Num_group+1;
        else
            loc2 = find(group(loc)==0);
            if loc2
                group(loc(loc2)) = maxval;
            end
            for com = minval:maxval
                if com > 0
                    loc2 = find(group(loc)==com);
                    if loc2
                        loc2 = find(group==com);
                        group(loc2) = maxval;
                    end
                end
            end
        end
    end
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
Num_group = max(group);


if Num_group>goal_num
    num = zeros(1,Num_group);
    for cnt =1: Num_group
        num(cnt) = length(find(group==cnt));
    end
    [snum IN] = sort(num,'descend');
    snum = snum(1:goal_num);
    locd = [];
    for cnt = 1:goal_num
        locd = [locd find(group==IN(cnt))];
    end
    for cnt = 1:Num_group
        if snum(goal_num)>num(cnt)
            locs = find(group==cnt); 
          %  locd = find(group~=cnt);
            if option == 'small'
                if 0
                    minval = min(min(W(locs,locd)));
                    pos = find(W(locs,locd)==minval);
                    pos1 = mod(pos(1),length(locs));
                    pos = (pos(1)-pos1)/length(locs)+1;
                else
                    [minval pos] = min(W(locs(1),locd));
                end
            else
                if 0
                    minval = max(max(W(locs,locd)));
                    pos = find(W(locs,locd)==minval);
                    pos1 = mod(pos(1),length(locs));
                    pos = (pos(1)-pos1)/length(locs)+1;
                else
                    [maxval pos] = max(W(locs(1),locd));
                end
            end
            group(locs) = group(locd(pos(1)))*ones(1,length(locs));
        end
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
    Num_group = max(group);
end


                

        
            