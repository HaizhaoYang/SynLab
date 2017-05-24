function [cluster, fq, amplt, num_freq] = freq_selection(T_f, is_real, eps, C, max_num, thre, res, pct, freq_range)
% Tracking an instantaneous frequency
%
% Get at most max_val components in the T_f distribution;
% Starting from the one with the largest T_f energy;
% Once one component has been identify, remove it from the T_f distribution
%
% By Haizhao Yang

[szm szn] = size(T_f);
loc = find(T_f<thre);
T_f(loc) = 0;
cluster = cell(max_num,1);
fq = [];
amplt = [];
Largest = max(max(T_f));
Largest = Largest(1);
for cnt = 1:max_num
    cls= [];
    max_val = max(max(T_f));
    max_val = max_val(1);
    if max_val < pct*Largest
        max_num = cnt -1;
        break;
    end
    [pos1,pos2] = find(T_f==max_val);
    pos1 = pos1(1); pos2 = pos2(1);
    R = round((pos1+freq_range(1)/res)*eps*C);
    if pos1-R>0
        st = pos1-R;
    else
        st  = 1;
    end
    if pos1+R<=szm
        ed = pos1+R;
    else
        ed = szm;
    end
    fq_range = st:ed;
    temp = T_f(fq_range,pos2);
    fq_temp = zeros(1,szn);
    amplt_temp = zeros(1,szn);
    fq_temp(pos2) = sum(temp.*fq_range')/sum(temp);
    amplt_temp(pos2) = sum(T_f(fq_range,pos2));
    loc = find(temp>0);
    cls = [cls fq_range(loc)+szm*(pos2-1)];
    T_f(fq_range,pos2) = 0;
    pos1 = round(fq_temp(pos2));
    if pos2~= 1
        cnt1 = pos1;
        for cnt2 = pos2-1:-1:1
            if cnt1-R>0
                st = cnt1-R;
            else
                st  = 1;
            end
            if cnt1+R<=szm
                ed = cnt1+R;
            else
                ed = szm;
            end
            fq_range = st:ed;
            temp = T_f(fq_range,cnt2);
            fq_temp(cnt2) = sum(temp.*fq_range')/sum(temp);
            amplt_temp(cnt2) = sum(T_f(fq_range,cnt2));
            loc = find(temp>0);
            cls = [cls fq_range(loc)+szm*(cnt2-1)];
            T_f(fq_range,cnt2) = 0;
            cnt1 = round(fq_temp(cnt2));
        end
    end
    if pos2~=szn
        cnt1 = pos1;
        for cnt2 = pos2+1:szn
            if cnt1-R>0
                st = cnt1-R;
            else
                st  = 1;
            end
            if cnt1+R<=szm
                ed = cnt1+R;
            else
                ed = szm;
            end
            fq_range = st:ed;
            temp = T_f(fq_range,cnt2);
            fq_temp(cnt2) = sum(temp.*fq_range')/sum(temp);
            amplt_temp(cnt2) = sum(T_f(fq_range,cnt2));
            loc = find(temp>0);
            cls = [cls fq_range(loc)+szm*(cnt2-1)];
            T_f(fq_range,cnt2) = 0;
            cnt1 = round(fq_temp(cnt2));
        end
    end
    fq = [fq; fq_temp];
    amplt = [amplt; amplt_temp];
    cluster{cnt} = cls;
end
num_freq = max_num;
if ~is_real
    fq = fq - szm/2;
end
fq = fq*res;
if freq_range(1)>0
    fq = fq + freq_range(1)-1;
end
cluster = cluster(1:num_freq,1);










