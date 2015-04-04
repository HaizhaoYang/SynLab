function [phase init] = phase_est(mode,pre_phase,amplt,is_real)
% Estimate phase functions 

[num_mode grid] = size(mode);
phase = zeros(num_mode,grid);
init = zeros(num_mode,1);
if is_real
    for cnt = 1:num_mode
        phase(cnt,:) = mono_arccos(mode(cnt,:)./amplt(cnt,:));
        init(cnt) = mean(phase(cnt,:)-pre_phase(cnt,:));
        phase(cnt,:) = pre_phase(cnt,:) +init(cnt);
    end
else
    phase = pre_phase;
end