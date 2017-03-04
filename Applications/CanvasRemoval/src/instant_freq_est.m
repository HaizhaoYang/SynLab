function [ins_fq spec] = instant_freq_est(fq,lowFq,is_mean)
num_fq = length(fq);
if num_fq > 1
    if lowFq >= fq(1)
        lowFq = round(fq(1)/4);
    end
    fq = sort(fq);
    max_trial = fq(1);
    dis = zeros(1,max_trial-lowFq);
    for cnt = lowFq+1:max_trial
        temp_ins_fq = fq(1)/cnt;
        temp = 0;
        for cnt2 = 2:num_fq
            int = fq(cnt2)/temp_ins_fq;
            temp = temp + (abs(int-round(int)));
        end
        dis(cnt-lowFq) = temp/(num_fq-1);
    end
    % dis
    [min_val min_pos] = min(dis);
    min_pos = min_pos + lowFq;
    spec = zeros(1,num_fq);
    for cnt = 1:num_fq
        spec(cnt) = round(fq(cnt)/min_pos);
    end
    if is_mean
        ins_fq = mean(fq/min_pos./spec')*min_pos;
    else
        ins_fq = fq(1);
    end
else
    ins_fq = fq;
    spec = 1;
end

