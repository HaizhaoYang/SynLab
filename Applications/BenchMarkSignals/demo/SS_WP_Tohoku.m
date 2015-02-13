clear all;
close all;


%% set up data
load data/tohoku.txt;
N = length(tohoku);
sampleRate = 20;
timeEnd = N/sampleRate;
upBound = round(sampleRate/2);
lowBound = 0;

x = [0:N-1]/N;
fff = tohoku(1:N);
fff = fff(:).';
fff = fff/max(abs(fff));

figure;plot(real(fff));

% parameters by default
R_high = round(N/10);       % range of interest is [R_low, R_high]
R_low = 0;
is_cos = 1;
is_unif = 1;
typeNUFFT = 1;
red = 64;
NG = round(N/2);    % number of subsampling points in space
epsl = 1e-10;         % threshold for SST
is_real = 1;        % 1: real signals, 0: complex signals

%% SS wave packet transform
if 1
    % parameters to be tuned
    res = 1;             % visualization resolution parameter in frequency
    rad = 1.2;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 1/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[0 upBound/20],ss_energy(1:end/4,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    colormap (1-gray); caxis([0,100]);
    saveas(pic,'results/Tohoku_SS_wp_1.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tohoku_SS_wp_1';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% SS wave packet transform
if 1
    % parameters to be tuned
    res = 1;             % visualization resolution parameter in frequency
    rad = 0.8;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 2/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[0 upBound/20],ss_energy(1:end/4,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    colormap (1-gray); caxis([0,100]);
    saveas(pic,'results/Tohoku_SS_wp_2.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tohoku_SS_wp_2';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% SS wavelet transform
if 1
    % parameters to be tuned
    res = 1;             % visualization resolution parameter in frequency
    rad = 0.5;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,1);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[0 upBound/20],ss_energy(1:end/4,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    colormap (1-gray); caxis([0,170]);
    saveas(pic,'results/Tohoku_SS_wl.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tohoku_SS_wl';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% SS wavelet transform
if 0
    % parameters to be tuned
    res = 1;             % visualization resolution parameter in frequency
    rad = 0.75;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[0 upBound/20],ss_energy(1:end/4,:));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;
    colormap (1-gray); caxis([0,0.05]);
    saveas(pic,'results/Tohoku_SS_wl.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Tohoku_SS_wl';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

