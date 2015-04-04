function phase = pre_phase_est(fq,dt)
% Integrate fq to get phase function

phase = zeros(size(fq));
for cnt = 2:size(fq,2)
    phase(:,cnt) = phase(:,cnt-1) + dt*fq(:,cnt);
end
    
