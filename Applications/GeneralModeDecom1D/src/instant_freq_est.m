function [ins_fq spec pre_ins_phase] = instant_freq_est(fq,dt,is_mean)
% Estmate the principle instantaneous frequency q(t) from a set of instantaneous
% freqencies { n_1 q(t), n_2 q(t), ... , n_k q(t)}
%
% By Haizhao Yang

[num_fq grid]= size(fq); 
if num_fq > 1
    fq = sortrows(fq,1);
    [min_val min_pos] = min(fq(1,:));
    max_trial = floor(min_val);
    dis = zeros(1,max_trial);
    for cnt = 1:max_trial
        temp_ins_fq = fq(1,:)/cnt;
        temp = 0;
        for cnt2 = 2:num_fq
            int = fq(cnt2,:)./temp_ins_fq;
            temp = temp + sum(abs(int-round(int)));
        end
        dis(cnt) = temp/length(grid)/(num_fq-1)/grid;
    end
   % dis
    [min_val min_pos] = min(dis)
    ins_fq = fq(1,:)/min_pos;
    spec = zeros(1,num_fq);
    for cnt = 1:num_fq
        spec(cnt) = round(mean(fq(cnt,:)./ins_fq));
    end
    if is_mean
        ins_fq = mean(fq./(spec'*ones(1,grid)));
    else
        ins_fq = fq(1,:);
    end
    pre_ins_phase = pre_phase_est(ins_fq,dt);
else
    ins_fq = fq;
    spec = 1;
    pre_ins_phase = pre_phase_est(ins_fq,dt);
end

    