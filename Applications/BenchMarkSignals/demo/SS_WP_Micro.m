clear all;
close all;

%% set up data
load data/microseismic.txt;
N = length(microseismic);

x = [0:N-1]/N;
fff = microseismic(1:N);
fff = fff(:).';
fff = fff/max(fff);

figure;plot(real(fff));

% parameters by default
epsl = 1e-10;         % threshold for SST
res = 0.5;             % visualization resolution parameter in frequency
NG = round(N);    % number of subsampling points in space
is_real = 1;        % 1: real signals, 0: complex signals       
R_low = 0;
is_cos = 1;
is_unif = 1;
typeNUFFT = 1;
red = 32;
sampleRate = 2000;
timeEnd = N/sampleRate;
upBound = round(sampleRate/2);
lowBound = 0;

%% SS wave packet transform
if 1
    % parameters to be tuned
    R_high = N/2;
    rad = 1.2;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 1/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic=figure;imagesc([0 timeEnd],[lowBound upBound],(ss_energy).^(1/3));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.7]);
    saveas(pic,'results/Micro_SS_wp_1.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_1';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % low frequency part
    R_high = round(R_high/2);
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic=figure;imagesc([0 timeEnd],[lowBound upBound/2],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.05]);
    saveas(pic,'results/Micro_SS_wp_1l.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_1l';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% SS wave packet transform
if 1
    % parameters to be tuned
    R_high = N/2;
    rad = 0.75;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 2/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[lowBound upBound],(ss_energy.^(1/3)));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.7]);
    saveas(pic,'results/Micro_SS_wp_2.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_2';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % low frequency part
    R_high = round(R_high/2);
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic=figure;imagesc([0 timeEnd],[lowBound upBound/2],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.05]);
    saveas(pic,'results/Micro_SS_wp_2l.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_2l';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end

%% SS wave packet transform
if 1
    % parameters to be tuned
    R_high = N/2;
    rad = 0.6;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 3/8;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[lowBound upBound],(ss_energy).^(1/3));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.7]);
    saveas(pic,'results/Micro_SS_wp_3.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_3';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % low frequency part
    R_high = round(R_high/2);
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic=figure;imagesc([0 timeEnd],[lowBound upBound/2],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.04]);
    saveas(pic,'results/Micro_SS_wp_3l.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wp_3l';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end


%% SS wavelet transform
if 1
    % parameters to be tuned
    R_high = N/2;
    rad = 0.4;          % rad in [0,2] to contral the size of supports of wave packets in the frequency domain
    t_sc = 1/2 + 1/2;        % t_sc in [0.5,1], scaling parameter, t_sc = 1 means wavelets, t_sc = 0.5 means wave atoms
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic = figure;imagesc([0 timeEnd],[lowBound upBound],(ss_energy).^(1/3));axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.7]);
    saveas(pic,'results/Micro_SS_wl.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wl';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
    
    % low frequency part
    R_high = round(R_high/2);
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,x,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
    
    % show results
    pic=figure;imagesc([0 timeEnd],[lowBound upBound/2],ss_energy);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;;
    colorbar off;colormap (1-gray);caxis([0 0.02]);
    saveas(pic,'results/Micro_SS_wl_l.fig');
    set(gca, 'FontSize', 18);
    b=get(gca);
    set(b.XLabel, 'FontSize', 18);set(b.YLabel, 'FontSize', 18);set(b.ZLabel, 'FontSize', 18);set(b.Title, 'FontSize', 18);
    str = 'results/Micro_SS_wl_l';
    print(gcf, '-depsc2', str);      command = sprintf('epstopdf %s.eps',str);      system(command);
end