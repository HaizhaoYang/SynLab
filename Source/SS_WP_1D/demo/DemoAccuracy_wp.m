% This code tests 1D (synchrosqueezed) fwd and inv wave packet transform with non-uniform data
% and arbitrary number of samples
%
% By Haizhao Yang

clear all;
close all;

%% Set up tests
test = cell(30,1); rslt = cell(30,1);
test{1,1} = struct('num',0,'is_real',0,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{2,1} = struct('num',0,'is_real',1,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{3,1} = struct('num',10,'is_real',0,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{4,1} = struct('num',10,'is_real',1,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{5,1} = struct('num',11,'is_real',0,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{6,1} = struct('num',11,'is_real',1,'is_unif',1,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{7,1} = struct('num',0,'is_real',0,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{8,1} = struct('num',0,'is_real',1,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{9,1} = struct('num',10,'is_real',0,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{10,1} = struct('num',10,'is_real',1,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{11,1} = struct('num',11,'is_real',0,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{12,1} = struct('num',11,'is_real',1,'is_unif',0,'typeNUFFT',1,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{13,1} = struct('num',0,'is_real',0,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{14,1} = struct('num',0,'is_real',1,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{15,1} = struct('num',10,'is_real',0,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{16,1} = struct('num',10,'is_real',1,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{17,1} = struct('num',11,'is_real',0,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{18,1} = struct('num',11,'is_real',1,'is_unif',0,'typeNUFFT',2,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{19,1} = struct('num',0,'is_real',0,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{20,1} = struct('num',0,'is_real',1,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{21,1} = struct('num',10,'is_real',0,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{22,1} = struct('num',10,'is_real',1,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{23,1} = struct('num',11,'is_real',0,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{24,1} = struct('num',11,'is_real',1,'is_unif',0,'typeNUFFT',3,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{25,1} = struct('num',0,'is_real',0,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{26,1} = struct('num',0,'is_real',1,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{27,1} = struct('num',10,'is_real',0,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{28,1} = struct('num',10,'is_real',1,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{29,1} = struct('num',11,'is_real',0,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);
test{30,1} = struct('num',11,'is_real',1,'is_unif',0,'typeNUFFT',4,'wpterr',[],'sstCerr',[],'sstTerr',[]);

%% Begin tests
for cnt = 1:6
    para = test{cnt,1};
    num = para.num;
    is_real = para.is_real;
    is_unif = para.is_unif;
    typeNUFFT = para.typeNUFFT;
    
    %set up data
    N = 1024-num;
     % assume the signal is defined on [0,1]
    xo = [0:N-1]/N;
    amp = 0.1;
    F1 = 60;
    F2 = 45;
    
    % xx is the nonlinear phase function of the first component
    xx = xo + amp*sin(2*pi*xo);
    f1u = 2*exp(2*pi*i*F1*2*xx);
    % yy is the nonlinear phase function of the second component
    yy = xo +amp*cos(2*pi*xo);
    f2u = exp(2*pi*i*F2*yy);
    % add noise if necessary
    NM = 0;
    ns = NM*randn(1,N);
    % obtain a superposition of two components
    ffu = f1u +f2u+ ns;
    
    if ~is_unif
        xo = xo + randn(1,N)/N/16;
        
        xx = xo + amp*sin(2*pi*xo);
        f1 = 2*exp(2*pi*i*F1*2*xx);
        yy = xo +amp*cos(2*pi*xo);
        f2 = exp(2*pi*i*F2*yy);
        
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
    % set up parameters, please see the comments in the corresponding 
    % source code for their definitions
    res = 2;
    NG = round(N/32);
    R_high = N/2;
    R_low = 0;
    rad = 1;
    is_cos = 1;
    t_sc = 1/2+1/4;
    epsl = 1e-4;
    red = 4;
    
    % perform the 1D forward wave packet transform
    coefCell = wp1_fwd(ffi,is_real,is_unif,typeNUFFT,xi,NG,R_high,R_low,rad,is_cos,t_sc,red);
    % coefCell is a cell structure storing the wave packet coefficients,
    % see the explanation in wp1_fwd for the arrangement of the cells
    
    % one could modify the content in coefCell for feature extraction,
    % signal processing, and then use the inverse transform below to see
    % the extracted signal or processed signal.
    
    % perform the 1D inverse wave packet transform
    ffr = wp1_inv(coefCell,N,is_real,R_high,R_low,rad,is_cos,t_sc);
    
    % visualize the results of the wave packet transform
    para.wpterr = [max(abs(ffr-ffc)) norm(ffc-ffr)/norm(ffc)];
    figure; 
    subplot(3,4,1);hold on;plot(real(ffr),'b');plot(real(ffc),'r');hold off;axis tight;
    subplot(3,4,2);hold on;plot(imag(ffr),'b');plot(imag(ffc),'r');hold off;axis tight;
    subplot(3,4,3);plot(real(ffr-ffc));axis tight;
    subplot(3,4,4);plot(imag(ffr-ffc));axis tight;
    
    % perform the 1D forward SS wave packet transform
    [ss_energy coefTensor InsFreq] = ss_wp1_fwd(ffi,is_real,is_unif,typeNUFFT,xi,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res);
    % ss_energy is a matrix storing the time-frequency representation by SS wave
    % packet transform
    % coefTensor is a cell structure storing the wave packet transform
    % coefficients
    % InsFreq is a cell structure storing the instantaneous frequency estimation
    % by the wave packet transform
    
    % there are two kinds of inverse SS wave packet transforms
    
    % the first kind is based on the features stored in a cell structure,
    % say the variable cluster below (1 for keeping and 0 for discarding); 
    % please see the input comments in ss_wp1_inv to see how to define the 
    % cluster variable.
    % in the example below, we do not modify the SS representation and
    % directly invert the result as follows.
    cluster = cell(1,1);  temp = zeros(size(ss_energy)); cluster{1} = find(temp==0);
    
    [ffr, amplt, num_mode] = ss_wp1_inv(cluster,coefTensor,InsFreq,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
    para.sstCerr = [max(abs(ffr-ffc)) norm(ffc-ffr)/norm(ffc)];
    subplot(3,4,5);hold on;plot(real(ffr),'b');plot(real(ffc),'r');hold off;axis tight;
    subplot(3,4,6);hold on;plot(imag(ffr),'b');plot(imag(ffc),'r');hold off;axis tight;
    subplot(3,4,7);plot(real(ffr-ffc));axis tight;
    subplot(3,4,8);plot(imag(ffr-ffc));axis tight;
    
    % the first kind is based on the features stored in a cell structure,
    % say the variable T below storing the SS representation; 
    % please see the input comments in
    % ss_wp1_inv to see how to define the T variable.
    % in the example below, we do not modify the SS representation and
    % directly invert the result as follows.
    T = cell(1,1); T{1} = ss_energy;
    [ffr, amplt, num_mode] = ss_wp1_invT(T,coefTensor,InsFreq,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
    para.sstTerr = [max(abs(ffr-ffc)) norm(ffc-ffr)/norm(ffc)];
    subplot(3,4,9);hold on;plot(real(ffr),'b');plot(real(ffc),'r');hold off;axis tight;
    subplot(3,4,10);hold on;plot(imag(ffr),'b');plot(imag(ffc),'r');hold off;axis tight;
    subplot(3,4,11);plot(real(ffr-ffc));axis tight;
    subplot(3,4,12);plot(imag(ffr-ffc));axis tight;
    rslt{cnt,1} = para;
    
    % if you would like to modify the SS representation and do the
    % transform, for example,
    T = cell(1,2);
    % the first component for the first half of the time domain
    T{1} = ss_energy;
    T{1}(:,end/2+1:end) = 0;
    % the second component for the second half of the time domain
    T{2} = ss_energy;
    T{2}(:,1:end/2) = 0;
    [ffr, amplt, num_mode] = ss_wp1_invT(T,coefTensor,InsFreq,is_real,N,R_high,R_low,rad,is_cos,t_sc,res);
    figure;subplot(1,2,1);plot(real(ffr(1,:))); subplot(1,2,2);plot(real(ffr(2,:)));
end
err = zeros(30,7);
for cnt = 1:6
    err(cnt,:) = [cnt rslt{cnt,1}.wpterr rslt{cnt,1}.sstCerr rslt{cnt,1}.sstTerr];
end
err(1:6,:)
%save 'Accuracy.mat' err rslt;


