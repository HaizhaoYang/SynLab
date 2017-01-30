clear all;
close all;

% Provided exact instantaneous frequencies and instantaneous amplitudes,
% this code shows the idea of diffeomorphism based spectral analysis.
%
% By Haizhao Yang

if(1)
    %set up data
    N = 8192;
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.005;
    is_plot = 0;
    
    num_group = 2;
    ins_freq = zeros(num_group,N);
    ins_amplt = zeros(num_group,N);
    ins_pre_phase = zeros(num_group,N);
    
    F1 = 30;
    xx = x + amp*sin(2*pi*x);
    ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*F1;
    f1 = zeros(1,N);
    am = 1+0.05*sin(4*pi*x);
    f1 = am.*gen_shape(F1*xx,2);
    
    F2 = 30;
    yy = x + amp*cos(2*pi*x);
    ins_freq(2,:) = (1-amp*2*pi*sin(2*pi*x))*F2;
    f2 = zeros(1,N);
    bm = 1+0.1*sin(2*pi*x);
    f2 = bm.*gen_shape(F2*yy,3);
    
    NM = 0;
    ns = NM*randn(1,N);
    fff = f1 + f2 + ns;
    
    
    ins_amplt(1,:) = am;
    ins_amplt(2,:) = bm;
    
    ins_pre_phase(1,:) = (xx-xx(1))*F1;
    ins_pre_phase(2,:) = (yy-yy(1))*F2;
    
    figure; hold on; plot(x,ins_freq(1,:),'r');plot(x,ins_freq(2,:),'b');title('Instantaneous frequencies');hold off;
end


if (1)
    %use smooth diffeomorphism and FFT to estimate spectrum of each general
    %wave shape
    %begin diffeomorphism
    tol = 1e-9;
    numPeaks = 200;
    
    
    signal_diff = zeros(num_group,N);
    spec_signal_diff = zeros(num_group,N);
    general_shape_comp = zeros(num_group,N);
    spec = cell(num_group,1);
    dist = zeros(num_group,N);
    uniform_sample = x;
    resid = fff;
    count = 0;
    
    if is_plot
        figure;
    end
    while norm(resid)>norm(fff)*tol & count < numPeaks
        count = count + 1;
        max_val = zeros(num_group,1);
        max_pos = zeros(num_group,1);
        for cnt = 1:num_group
            nonuniform_sample = ins_pre_phase(cnt,:)/mean(ins_freq(cnt,:));
            signal_diff(cnt,:) = spline(nonuniform_sample,resid./ins_amplt(cnt,:),uniform_sample);
            %[T_f coef kk] = ss_wp_1D(signal_diff(cnt,:),1,freq_range,eps,res,NG,szchar);
            %figure;imagesc(T_f);head = sprintf('%d smooth deffeomorphism with pre_phase',cnt);title(head);
            %FFT
            spec_signal_diff(cnt,:) = ifftshift(fft(fftshift(signal_diff(cnt,:))))/length(signal_diff(cnt,:));
            if is_plot
                if cnt ==1
                    subplot(1,2,1);
                else
                    subplot(1,2,2);
                end
                plot(-N/2:(N-1)/2,abs(spec_signal_diff(cnt,:)));
            end
            [temp_val temp_pos] = max(abs(spec_signal_diff(cnt,1:N/2)));
            max_val(cnt,1) = temp_val(1);
            max_pos(cnt,1) = abs(temp_pos(1)-N/2-1);
        end
        if is_plot
            pause;
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
    if 0
        for cnt = 1:num_group
            temp = dist(cnt,N/2+2:end);
            loc = find(temp>0);
            spec{cnt,1} = unique(round(loc/loc(1)));
            spec{cnt,1}
        end
    end
    err1 = abs(real(f1-general_shape_comp(1,:)));
    err2 = abs(real(f2-general_shape_comp(2,:)));
    if is_plot
        if num_group >1
            figure;
            subplot(1,2,1);plot(dist(1,:));title('Spectrum of 1st general wave shape component');
            subplot(1,2,2);plot(dist(2,:));title('Spectrum of 2nd general wave shape component');
            figure;
            subplot(1,2,1);plot(err1);title('error of 1st general component');
            subplot(1,2,2);plot(err2);title('error of 2nd general component');
            figure;
            subplot(1,2,1); hold on; plot(f1(1:100),'b'); plot(general_shape_comp(1,1:100),'r'); hold off;
            subplot(1,2,2); hold on; plot(f2(1:100),'b'); plot(general_shape_comp(2,1:100),'r'); hold off;
        else
            figure;plot(general_shape_comp(1,1:100));
        end
    end
    pic = figure; plot(x(1:2000),f1(1:2000)+f2(1:2000));set(gca,'DataAspectRatio',[1 180 1]);xlabel('Time (Second)');title('A superposition of components');
    
    [val pos] = max(err1);
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    pic = figure; hold on; plot(x(posVec),general_shape_comp(1,posVec),'r'); plot(x(posVec),f1(posVec),'b');axis square;title('1st mode');
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    pic = figure; hold on; plot(x(posVec),general_shape_comp(2,posVec),'r'); plot(x(posVec),f2(posVec),'b');axis square;title('2nd mode');
end


