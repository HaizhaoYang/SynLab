function rsl = sharpen(binaryIndex,numSmooth,len,var)
if nargin<2, numSmooth = 20; end;
if nargin<3, len = 5; end;
if nargin<4, var = 1; end;
temp = smoothImage(binaryIndex,5,1);
pos = find(binaryIndex==1);
for cnt = 1:numSmooth
    temp = smoothImage(temp,len,var);
    temp(pos) = 1;
end
tempp = temp;
[sm sn]=size(temp);
direction = 1;
for cnt1 = 1:sm
    for cnt2 = 1:sn
        if cnt2 < sn
            if direction == 1
                if temp(cnt1,cnt2) >= temp(cnt1,cnt2+1)
                    temp(cnt1,cnt2) = 3;
                else
                    temp(cnt1,cnt2) = 0;
                    direction = -direction;
                end
            else
                if temp(cnt1,cnt2) <= temp(cnt1,cnt2+1)
                    temp(cnt1,cnt2) = 3;
                else
                    direction = -direction;
                end
            end
        else
            if direction == 1
                if temp(cnt1,cnt2) >= temp(cnt1,1)
                    temp(cnt1,cnt2) = 3;
                else
                    temp(cnt1,cnt2) = 0;
                    direction = -direction;
                end
            else
                if temp(cnt1,cnt2) <= temp(cnt1,1)
                    temp(cnt1,cnt2) = 3;
                else
                    direction = -direction;
                end
            end
        end
    end
end
pos = find(temp>0);
temp(pos) = 3;

hon = temp;
temp = tempp;
direction = -1;
for cnt2 = 1:sn
    for cnt1 = 1:sm
        if cnt1 < sm
            if direction == 1
                if temp(cnt1,cnt2) >= temp(cnt1+1,cnt2)
                    temp(cnt1,cnt2) = 3;
                else
                    temp(cnt1,cnt2) = 0;
                    direction = -direction;
                end
            else
                if temp(cnt1,cnt2) <= temp(cnt1+1,cnt2)
                    temp(cnt1,cnt2) = 3;
                else
                    direction = -direction;
                end
            end
        else
            if direction == 1
                if temp(cnt1,cnt2) >= temp(1,cnt2)
                    temp(cnt1,cnt2) = 3;
                else
                    temp(cnt1,cnt2) = 0;
                    direction = -direction;
                end
            else
                if temp(cnt1,cnt2) <= temp(1,cnt2)
                    temp(cnt1,cnt2) = 3;
                else
                    direction = -direction;
                end
            end
        end
    end
end
pos = find(temp>0);
temp(pos) = 3;

ver = temp;

rsl = min(hon,ver);



