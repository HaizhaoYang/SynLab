% This code tests the time complexity of 1D (synchrosqueezed) fwd and inv
% wave packet transform with arbitrary number of samples
%
% By Haizhao Yang

clear all;
close all;

%% Set up tests
test = cell(1,1); rslt = cell(1,1);
test{1,1} = struct('num',0,'is_real',0,'is_unif',1,'typeNUFFT',1,'wpttime',[],'sstime',[],'isstime1',[],'isstime2',[]);

%% Begin tests
for cnt = 1:1
    for pow = 20
        test{cnt,1} = test{cnt,1};
        num = test{cnt,1}.num;
        is_real = test{cnt,1}.is_real;
        is_unif = test{cnt,1}.is_unif;
        typeNUFFT = test{cnt,1}.typeNUFFT;
        
        N = 2^pow-num;
        xo = [0:N-1]/N;
        amp = 0.1;
        F1 = 20;
        F2 = 15;
        
        xx = xo + amp*sin(2*pi*xo);
        f1u = [0:N-1]/N;%2*exp(2*pi*i*F1*2*xx);
        yy = xo +amp*cos(2*pi*xo);
        f2u = [0:N-1]/N;%exp(2*pi*i*F2*yy);
        
        NM = 0;
        ns = NM*randn(1,N);
        ffu = f1u +f2u+ ns;
        
        if ~is_unif
            xo = xo + randn(1,N)/N/16;
            
            xx = xo + amp*sin(2*pi*xo);
            f1 = [0:N-1]/N;%2*exp(2*pi*i*F1*2*xx);
            yy = xo +amp*cos(2*pi*xo);
            f2 = [0:N-1]/N;%exp(2*pi*i*F2*yy);
            
            NM = 0;
            ns = NM*randn(1,N);
            fff = f1 +f2+ ns;
        else
            fff = ffu;
        end
        
        if is_real
            ffi = real(fff);
            ffc = real(ffu);
        else
            ffi = fff;
            ffc = ffu;
        end
        xi = xo;
        res = 2;
        NG = 1;%round(N/32);
        R_high = N/2;
        R_low = 0;
        rad = 1;
        is_cos = 1;
        t_sc = 1/2;
        epsl = 1e-4;
        red = 1;
        
        tic;
        coefCell = wp1_fwd(ffi,is_real,is_unif,typeNUFFT,xi,NG,R_high,R_low,rad,is_cos,t_sc,red);
        ffr = wp1_inv(coefCell,N,is_real,R_high,R_low,rad,is_cos,t_sc);
        test{cnt,1}.wpttime = [test{cnt,1}.wpttime toc];
%         tic;
%         [ss_energy coefTensor InsFreq] = ss_wp1_fwd(ffi,is_real,is_unif,typeNUFFT,xi,NG,R_high,R_low,rad,is_cos,t_sc,epsl,res);
%         test{cnt,1}.sstime = [test{cnt,1}.sstime toc];
% %         cluster = cell(1,1);  temp = zeros(size(ss_energy)); cluster{1} = find(temp==0);
% %         tic;
% %         [ffr, amplt, num_mode] = ss_wp1_inv(cluster,coefTensor,InsFreq,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
% %         test{cnt,1}.isstime1 = [test{cnt,1}.isstime1 toc];
%         T = cell(1,1); T{1} = ss_energy;
%         tic;
%         [ffr, amplt, num_mode] = ss_wp1_invT(T,coefTensor,InsFreq,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
%         test{cnt,1}.isstime2 = [test{cnt,1}.isstime2 toc];
    end
end


