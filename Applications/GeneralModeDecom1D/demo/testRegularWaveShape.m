clear all;
close all;

% This code solves the general mode decomposition problem under the weak
% well-separation and well-different condition in the paper "Synchrosqueezed
% wave packet transforms and diffeomorphism based spectral analysis for 1D 
% general mode decompositions", H. Yang, Applied and Computational Harmonic
% Analysis, 2014.
%
% Test general wave shape functions 
%
% By Haizhao Yang

if(1)
    %set up data
    N = 8192;
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.01;
    F1 = 60;
    F2 = 90;
    xx = x + amp*sin(2*pi*x);
    f1 = zeros(1,N);
    am = 1+0.05*sin(4*pi*x);
    f1 = gen_shape2(F1*xx,4);%4
    yy = x + amp*cos(2*pi*x);
    %yy = x + amp*sin(4*pi*x);
    f2 = zeros(1,N);
    bm = 1+0.1*sin(2*pi*x);
    f2 = gen_shape(F2*yy,1);%5
    
    NM = 0;% is fine
    ns = NM*randn(1,N);
    fff = f1 + f2 + ns;
end

if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-3;
    res =2;
    freq_range = [0 N/2];
    NG = N/4;
    dt = 1/NG;
    t_sc = 0.75;
    [T_f coef kk] = ss_wp1_fwd(fff,1,1,1,x,NG,freq_range(2),freq_range(1),1,1,t_sc,1,eps,res);
    loc = find(T_f<NM^2);
    T_f(loc) = 0;
    pic = figure;imagesc([0 1],[1 600],T_f(1:300,:));colormap (1-gray);axis square;xlabel('Time (Second)');ylabel('Frequency (Hz)');axis xy;title('SSWPT');
end


if (1)
    %synchrosqueezed wave packet transform
    eps = 1e-3;
    res =1;
    freq_range = [0 N/2];
    NG = N;
    dt = 1/NG;
    t_sc = 0.75
    [T_f coef kk] = ss_wp1_fwd(fff,1,1,1,x,NG,freq_range(2),freq_range(1),1,1,t_sc,1,eps,res);
    loc = find(T_f<NM^2);
    T_f(loc) = 0;
end

if (1)
    %extract isolated component and select instantaneous freqencies
    thre = 0;
    C = 5;
    max_num = 3;
    pct = 0.01;
    Bound = 250;
    T_f_select = zeros(size(T_f));
    T_f_select(1:Bound,:) = T_f(1:Bound,:);
    [cluster, freq, amplt, num_freq] = freq_selection(T_f_select, 1, eps*20, C, max_num, thre, res, pct, freq_range);
    comp_select = zeros(num_freq,N);
    pre_amplt = zeros(num_freq,N);
    for cnt = 1:num_freq
        freq(cnt,:) = filt_hard_threshold(freq(cnt,:),0.01,1);
        T_temp = cell(1,1);
        T_temp{1} = zeros(size(T_f));
        T_temp{1}(cluster{cnt,1}) = T_f(cluster{cnt,1});
        [comp_select(cnt,:), num] = ss_wp1_invT(T_temp, coef, kk, 1, N, freq_range(2), freq_range(1), 1, 1, t_sc, res); 
        pre_amplt(cnt,:) = amplt_est(1,comp_select(cnt,:));  
        pre_amplt(cnt,:) = filt_hard_threshold(pre_amplt(cnt,:),0.01,1);
    end
    %get exact frequency and phase
    phase = zeros(num_freq,N);
    for cnt = 1:num_freq
        pre_phase(cnt,:) = pre_phase_est(freq(cnt,:),dt);
    end
    freq = filterFreq(freq,0);
    %----------------------------------------------------------------------
    %use linear regression to quantify correlation
    %use spectral clustering to separate
    sigma_slop = 0.01;
    sigma_resid = 1.5;
    Amat_slope = zeros(num_freq,num_freq);
    Amat_resid = zeros(num_freq,num_freq);
    for cnt1 = 1:num_freq-1
        for cnt2 = cnt1+1:num_freq
            [dis resid] = distance_linear_fit(freq(cnt1,:),freq(cnt2,:));
            Amat_slope(cnt1,cnt2) = gauss(dis,sigma_slop);
            Amat_resid(cnt1,cnt2) = gauss(resid,sigma_resid);
        end
    end
    Amat_slope = max(Amat_slope,Amat_slope') + eye(num_freq);
    Amat_resid = max(Amat_resid,Amat_resid') + eye(num_freq);
    if num_freq > 3
        [group_resid , svals, Lap , num_group_resid] = SpectralClustering_est(Amat_resid,3);
    else if num_freq == 1
            group_resid = [1];
            num_group_resid = 1;
        else
        if Amat_resid(1,2) > 0.5
            num_group_resid = 1;
            group_resid = [1 1];
        else
            num_group_resid = 2;
            group_resid = [1 2];
        end
        end
    end
    group = group_resid;
    num_group = num_group_resid;
    %----------------------------------------------------------------------
    %search for instantaneous frequency, ins_pre_amplt,ins_pre_phase
    ins_freq = zeros(num_group,N);
    ins_amplt = zeros(num_group,N);
    ins_pre_phase = zeros(num_group,N);
    pre_spec = cell(1,num_group);
    for cnt = 1:num_group
        loc = find(group == cnt);
        temp_fq = freq(loc,:);
        [ins_freq(cnt,:) pre_spec{cnt}] = instant_freq_est(temp_fq,dt);
        [max_val max_pos] = max(sum(pre_amplt(loc,:)'));
        ins_amplt(cnt,:) = pre_amplt(loc(max_pos),:);
        loc2 = find(pre_spec{cnt}==1);
        ins_pre_phase(cnt,:) = pre_phase(loc(loc2),:);
    end
    %show result
    if 0
        for cnt = 1:num_group
            head = sprintf('%d instantaneous frequency',cnt);
            figure; plot(ins_freq(cnt,:));title(head);
            head = sprintf('%d instantaneous amplitude',cnt);
            figure; plot(ins_amplt(cnt,:));title(head);   
        end
    end
end


if (1)
    %use smooth diffeomorphism and FFT to estimate spectrum of each general
    %wave shape
    %begin diffeomorphism
    %use pre_phase
    signal_diff = zeros(num_group,N);
    spec_signal_diff = zeros(num_group,N);
    general_shape_comp = zeros(num_group,N);
    spec = cell(num_group,1);
    dist = zeros(num_group,N);
    uniform_sample = x;
    resid = fff;
    count = 0;
    while norm(resid)>norm(fff)*0.0000001 & count < 50
        count = count + 1;
        max_val = zeros(num_group,1);
        max_pos = zeros(num_group,1);
        for cnt = 1:num_group
            nonuniform_sample = ins_pre_phase(cnt,:)/mean(ins_freq(cnt,:));
            signal_diff(cnt,:) = spline(nonuniform_sample,resid./ins_amplt(cnt,:),uniform_sample);
            spec_signal_diff(cnt,:) = ifftshift(fft(fftshift(signal_diff(cnt,:))))/length(signal_diff(cnt,:));
            [temp_val temp_pos] = max(abs(spec_signal_diff(cnt,1:N/2)));
            max_val(cnt,1) = temp_val(1);
            max_pos(cnt,1) = abs(temp_pos(1)-N/2-1);
        end
        [max_val2 max_pos2] = max(max_val);
        max_val2 = max_val2(1);
        max_pos2 = max_pos2(1);
        dist(max_pos2,max_pos(max_pos2,1)+N/2+1) = 1;
        dist(max_pos2,N/2+1-max_pos(max_pos2,1)) = 1;
        %estimate component for given pre_phase and ins_amplt
        %max_pos(max_pos2,1)*pre_phase(max_pos2,:)/mean(freq(max_pos2,:))
        %The general wave shape is real, so comp corresponds to s(n) and
        %s(-n).
        [comp amplt_cos amplt_sin] = opt_amplt(resid,ins_amplt(max_pos2,:),max_pos(max_pos2,1)*ins_pre_phase(max_pos2,:)/mean(ins_freq(max_pos2,:)));
        general_shape_comp(max_pos2,:) = general_shape_comp(max_pos2,:) + comp;
        resid = resid - comp;
    end
    %estimate spectrum from dist
    for cnt = 1:num_group
        temp = dist(cnt,N/2+2:end);
        loc = find(temp>0);
        spec{cnt,1} = unique(round(loc/loc(1)));
        spec{cnt,1}
    end
    pic = figure; plot(x(1:2000),f1(1:2000)+f2(1:2000));set(gca,'DataAspectRatio',[1 180 1]);xlabel('Time (Second)');title('A superposition of general modes');
    pic = figure; hold on; plot(x,ins_freq(1,:),'r'); plot(x,(1+amp*2*pi*cos(2*pi*x))*F1,'b');
    plot(x,ins_freq(2,:),'r'); plot(x,(1-amp*2*pi*sin(2*pi*x))*F2,'b');axis square;
    xlabel('Time (Second)');ylabel('Frequency (Hz)');title('Instantaneous frequencies');
    pic = figure; hold on; plot(x(1:300),general_shape_comp(2,1:300),'r'); plot(x(1:300),f1(1:300),'b');axis square;title('1st mode');
    pic = figure; hold on; plot(x(1:300),general_shape_comp(1,1:300),'r'); plot(x(1:300),f2(1:300),'b');axis square;title('2nd mode');
    A = [fff(1:1500);f1(1:1500);general_shape_comp(2,1:1500);f2(1:1500);general_shape_comp(1,1:1500)];
    pic = plot_table(A,[-5 5]);
end
