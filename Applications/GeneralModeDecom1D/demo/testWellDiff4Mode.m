clear all;
close all;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provided exact instantaneous frequencies and instantaneous amplitudes,
% this code tests the idea of diffeomorphism based spectral analysis (DSA).
%
% Larger 'dif' or 'FF' leads to better 'well-difference' of instantaneous frequencies
% and hence the DSA method works better.
%
% Smaller 'dif' or 'FF' violates the 'well-difference' condition and the
% DSA method fails.
%
% By Haizhao Yang

%%
for dif = 1:-0.2:0
    %set up data
    N = 8192;
    x = [0:N-1]/N;
    fff = zeros(1,N);
    amp = 0.005;
    is_plot = 0;
    FF = 160;
    
    num_group = 4;
    ins_freq = zeros(num_group,N);
    ins_amplt = zeros(num_group,N);
    ins_pre_phase = zeros(num_group,N);
    
    xx = x + amp*sin(2*pi*x);
    ins_freq(1,:) = (1+amp*2*pi*cos(2*pi*x))*FF;
    f1 = zeros(1,N);
    ins_amplt(1,:) = 1+0.05*sin(4*pi*x);
    f1 = ins_amplt(1,:).*gen_shape2(FF*xx,3);
    
    yy = x + amp*sin(2*pi*x+dif);
    ins_freq(2,:) = (1+amp*2*pi*cos(2*pi*x+dif))*FF;
    f2 = zeros(1,N);
    ins_amplt(2,:) = 1+0.1*sin(2*pi*x);
    f2 = ins_amplt(2,:).*gen_shape(FF*yy,2);
    
    zz = x + amp*sin(2*pi*x+dif*2);
    ins_freq(3,:) = (1+amp*2*pi*cos(2*pi*x+dif*2))*FF;
    f3 = zeros(1,N);
    ins_amplt(3,:) = 1+0.2*sin(2*pi*x);
    f3 = ins_amplt(3,:).*gen_shape(FF*zz,1);
    
    ww = x + amp*sin(2*pi*x+dif*3);
    ins_freq(4,:) = (1+amp*2*pi*cos(2*pi*x+dif*3))*FF;
    f4 = zeros(1,N);
    ins_amplt(4,:) = 1+0.2*sin(2*pi*x);
    f4 = ins_amplt(4,:).*gen_shape(FF*ww,3);
    
    NM = 0;
    ns = NM*randn(1,N);
    fff = f1 + f2 + f3 + f4 + ns;
    
    ins_pre_phase(1,:) = (xx-xx(1))*FF;
    ins_pre_phase(2,:) = (yy-yy(1))*FF;
    ins_pre_phase(3,:) = (zz-zz(1))*FF;
    ins_pre_phase(4,:) = (ww-ww(1))*FF;
    
    %use smooth diffeomorphism and FFT to estimate spectrum of each general
    %wave shape
    %set up parameters
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
            spec_signal_diff(cnt,:) = ifftshift(fft(fftshift(signal_diff(cnt,:))))/length(signal_diff(cnt,:));
            if is_plot
                if cnt == 1
                    subplot(1,4,1);
                else if cnt == 2
                        subplot(1,4,2);
                    else if cnt == 3
                            subplot(1,4,3);
                        else
                            supblot(1,4,4);
                        end
                    end
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
    err3 = abs(real(f3-general_shape_comp(3,:)));
    err4 = abs(real(f4-general_shape_comp(4,:)));
    
    % show results
    figure; subplot(3,2,1);hold on; plot(x,ins_freq(1,:),'r');plot(x,ins_freq(2,:),'b');plot(x,ins_freq(3,:),'g');plot(x,ins_freq(4,:),'k');title('Instantaneous frequencies');hold off;
    subplot(3,2,2); plot(x(1:2000),fff(1:2000));set(gca,'DataAspectRatio',[1 180 1]);xlabel('Time (Second)');title('A superposition of components');
    
    [val pos] = max(err1);
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    subplot(3,2,3); hold on; plot(x(posVec),general_shape_comp(1,posVec),'r'); plot(x(posVec),f1(posVec),'b');axis square;title('1st mode');
    [val pos] = max(err2);
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    subplot(3,2,4); hold on; plot(x(posVec),general_shape_comp(2,posVec),'r'); plot(x(posVec),f2(posVec),'b');axis square;title('2nd mode');
    [val pos] = max(err3);
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    subplot(3,2,5); hold on; plot(x(posVec),general_shape_comp(3,posVec),'r'); plot(x(posVec),f3(posVec),'b');axis square;title('3rd mode');
    [val pos] = max(err4);
    posVec = pos-100:pos+100;
    posVec = posVec - min(posVec(1)-1,0);
    posVec = posVec - max(posVec(end)-N,0);
    subplot(3,2,6); hold on; plot(x(posVec),general_shape_comp(4,posVec),'r'); plot(x(posVec),f4(posVec),'b');axis square;title('4th mode');
    
    if is_plot
        figure;
        subplot(2,4,1);plot(dist(1,:));title('Spectrum of 1st general wave shape component');
        subplot(2,4,2);plot(dist(2,:));title('Spectrum of 2nd general wave shape component');
        subplot(2,4,3);plot(dist(3,:));title('Spectrum of 3rd general wave shape component');
        subplot(2,4,4);plot(dist(4,:));title('Spectrum of 4th general wave shape component');
        
        subplot(2,4,5);plot(err1);title('error of 1st general component');
        subplot(2,4,6);plot(err2);title('error of 2nd general component');
        subplot(2,4,7);plot(err3);title('error of 3rd general component');
        subplot(2,4,8);plot(err4);title('error of 4th general component');
    end
end

