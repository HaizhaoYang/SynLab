clear all;
close all;


%% set up data
load data/tremor.txt;
N_truncation = round(length(tremor)*0.9);
N = length(tremor);

x = [0:N-1]/N;
fff = tremor(1:N);
fff = fff(:).';
fff(N_truncation:end) = 0;
%fff = fff/max(fff);

figure;plot(real(fff));

% parameters by default
epsl = 1e-10;         % threshold for SST
R_high = 36000;
R_low = 0;
res = 30;%30;             % visualization resolution parameter in frequency
is_cos = 1;
is_real = 1;        % 1: real signals, 0: complex signals
is_unif = 1;
typeNUFFT = 1;
red = 100;
upBound = round(R_high/2/res);
lowBound = max(1,round(R_low/2/res));

%% synchrosqueezed wave packet transform
if 1
    % parameters to be tuned
    NG = round(N/120);    % number of subsampling points in space
    rad = 0.7;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 1/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    
    % SST
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    [m,n] = size(ss_energy);
    temp = ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    %pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    %colorbar off; caxis([0 0.01]);colormap (1-gray);
    
    % better visualization in the high frequency part
    ss_energy = (diag(((1:m).^2)))*ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off; caxis([0 0.02]);colormap (1-gray);
    saveas(pic,'results/Tremor_SS_wp_1.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tremor_SS_wp_1';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
end

%% synchrosqueezed wave packet transform
if 1
    % parameters to be tuned
    NG = round(N/30);
    rad = 0.4;
    t_sc = 1/2 + 2/8;
    
    % SST
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    [m,n] = size(ss_energy);
    temp = ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    %pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    %colorbar off; caxis([0 0.01]);colormap (1-gray);
    
    % better visualization in the high frequency part
    ss_energy = (diag(((1:m).^2)))*ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off; caxis([0 0.01]);colormap (1-gray);
    saveas(pic,'results/Tremor_SS_wp_2.fig');
    
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tremor_SS_wp_2';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% synchrosqueezed wavelet transform
if 1
    % parameters to be tuned
    NG = round(N/20);
    rad = 0.05;
    t_sc = 1;
    
    % SST
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    [m,n] = size(ss_energy);
    temp = ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    %pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    %colorbar off; caxis([0 0.01]);colormap (1-gray);
    
    % better visualization in the high frequency part
    ss_energy = (diag(((1:m).^2)))*ss_energy;
    ss_energy = ss_energy/max(ss_energy(:));
    
    % show results
    pic = figure;imagesc([0 600],[2*res*lowBound 2*res*upBound]/1200,ss_energy(lowBound:upBound,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off; caxis([0 0.005]);colormap (1-gray);
    saveas(pic,'results/Tremor_SS_wl.fig');
    
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tremor_SS_wl';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end
