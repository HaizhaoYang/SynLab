function [BD agl] = initialGuessOneSector(fff,subsampleRate,plotFIG,parameters)
%function [R agl TTEng_1st TTEng_2nd W_sec ss_energy] = initialGuess(fff,subsampleRate,plotFIG,parameters)
% This code uses ss curvelet transforms to detect different oscillation directions.
%% Main steps:
% 1) Detect useful frequency band;
% 2) Apply the synchrosqueezed transform SS_ct_polar() to obtain the ss energy
% distribution ss_energy
% 3) Apply bump detection algorithm Local_elastic() to estimate local wave
% vectors
%% Remark:
%1) At a spacial grid (x1,x2), the location of the largest energy bump in the i^th sector is stored in a polar
%coordinate with R(x1,x2,i) and agl(x1,x2,i);
%2) The energy of this largest energy bump is TTEng_1st(x1,x2,i);
%3) TTEng_2nd(x1,x2,i) is the energy of the 2nd largest energy bump;
%4) W_sec(x1,x2,i) = TTEng_1st(x1,x2,i) / total energy in the i^th
%sector at location (x1,x2);
%5) The variables ending with '_full' are full size data for the
%original crystal image; other variables are down sampled data and
%their sizes are controled by the variable 'subsampleRate';
%% set up parameters
if nargin<4
    parameters = struct('extendRatio',0.1,'threRatio',0.1,'epsFreqBand',1e-2,'numWave',3,'numGridPerSector',30,'red',2,'rad',2,'s_sc',0.8,'t_sc',0.8,'epsSST',1e-2);
end
if nargin<3, plotFIG = false; end;
if nargin<2, subsampleRate = 2; end;
%% normalize the image
fff = fff/max(max(abs(fff)));

% %% estimate Frequency band
% tic;
extendRatio = parameters.extendRatio;%influence the estimate around the boundaries.
threRatio = parameters.threRatio;
epsl = parameters.epsFreqBand;

if 1
    [m n] = size(fff);
    temp = zeros(n,m);
    for cnt1 = 1:n
        for cnt2 = 1:m
            temp(cnt1,cnt2) = fff(m+1-cnt2,cnt1);
        end
    end
    fff = temp;
end

% [m n] = size(fff);
% mm = min(m,n);
% fff = fff(1:mm-1,1:mm-1);
[m n] = size(fff);
mx = m;
nx = n;
subsampleRate = 4;
dm = round(m*0.2/subsampleRate)*subsampleRate;
dn = round(n*0.2/subsampleRate)*subsampleRate;
fff = [zeros(dm,dn*2+n);zeros(m,dn) fff zeros(m,dn);zeros(dm,dn*2+n)];


[m n] = size(fff);
mmx = m; nnx = n;
if m ~= n
    mm = max(m,n);
    if m > n
        type = 0;
        X = [fff zeros(m,m-n)];
    else
        type = 1;
        X = [fff; zeros(n-m,n)];
    end
    idx = 1:mm;
    X = X(idx,idx);
else 
    type = 0;
    X = fff;
    mm = m;
end
Y=fftshift(fft2(ifftshift(X)));
idx = -(mm+mod(mm+1,2)-1)/2:(mm-mod(mm+1,2)-1)/2;
[X1 X2] = ndgrid(idx,idx);
R = round(sqrt(X1.^2+X2.^2));
R = R(:);
Y = Y(:);
thre = max(Y)*epsl;
pos = find(Y>thre);
L = max(R);
spec = zeros(1,L+1);
for cnt = 1:length(pos)
    spec(R(pos(cnt))+1) = spec(R(pos(cnt))+1) + Y(pos(cnt))*Y(pos(cnt))';
end
spec(1) = 0;
[val pos] = max(spec);
thre = val*threRatio;
ed = pos+1;
cnt = 0;
while cnt == 0 & ed < L
    if spec(ed)>=spec(ed-1) & spec(ed) < thre
        cnt = 1;
    end
    ed = ed + 1;
end
ed = min(L,ceil((ed-1)*(1+extendRatio)));
st = pos -1;
cnt = 0;
while cnt == 0 & st > 1
    if spec(st)>=spec(st+1) & spec(st) < thre
        cnt = 1;
    end
    st = st - 1;
end
st = max(1,floor((st+1)*(1-extendRatio)));
if st <= 2
    spec(1:ed) = 0;
    [val pos] = max(spec);
    thre = val*threRatio;
    ed = pos+1;
    cnt = 0;
    while cnt == 0 & ed < L
        if spec(ed)>=spec(ed-1) & spec(ed) < thre
            cnt = 1;
        end
        ed = ed + 1;
    end
    ed = min(L,ceil((ed-1)*(1+extendRatio)));
    st = pos -1;
    cnt = 0;
    while cnt == 0 & st > 1
        if spec(st)>=spec(st+1) & spec(st) < thre
            cnt = 1;
        end
        st = st - 1;
    end
    st = max(1,floor((st+1)*(1-extendRatio)));
end
if plotFIG
    figure;hold on; plot(spec); plot(st,spec(st),'ro'); plot(ed,spec(ed),'ro');hold off;
    title('Fourier energy distribution and identified energy band');
end
%% Apply synchrosqueezed transforms
R_low = st; R_high = ed;
%set up the size of samples in [-R_high,R_high]*[-R_high,R_high] in the
%phase domain
num_wave = parameters.numWave;%determined by the reference grain
NB = [1,num_wave*parameters.numGridPerSector];%[(R_high - R_low),num_wave*parameters.numGridPerSector];
%set up parameters for ss transform
red = parameters.red;
rad = parameters.rad;
epsl = parameters.epsSST;
t_sc = parameters.t_sc;
s_sc = parameters.s_sc;
num_direction = 1;
is_real = parameters.is_real;

m = mx;
n = nx;
R_low = st; R_high = ed;
SPg = round(size(X)/subsampleRate);

ss_energy = SS_ct_polar_v1(num_direction,X,SPg,NB,rad,is_real,R_low,R_high,epsl,red(1),t_sc, s_sc);
[agl TTEng_1st TTEng_2nd] = LocSmooth(sum(ss_energy),num_wave);
BD = 1./sqrt(TTEng_1st-TTEng_2nd+1);

if type == 0  
    BD = BD(:,1:round(nnx/subsampleRate));
else
    BD = BD(1:round(mmx/subsampleRate),:);
end
BD = BD(dm/subsampleRate+1:round((m+dm)/subsampleRate),dn/subsampleRate+1:round((n+dn)/subsampleRate));
if 1
    [tempm tempn] = size(BD);
    temp = zeros(tempn,tempm);
    for cnt1 = 1:tempn
        for cnt2 = 1:tempm
            temp(cnt1,cnt2) = BD(cnt2,tempn+1-cnt1);
        end
    end
    BD = temp;
end

if type == 0  
    agl = agl(:,1:round(nnx/subsampleRate));
else
    agl = agl(1:round(mmx/subsampleRate),:);
end
agl = agl(dm/subsampleRate+1:round((m+dm)/subsampleRate),dn/subsampleRate+1:round((n+dn)/subsampleRate));
if 1
    [tempm tempn] = size(agl);
    temp = zeros(tempn,tempm);
    for cnt1 = 1:tempn
        for cnt2 = 1:tempm
            temp(cnt1,cnt2) = agl(cnt2,tempn+1-cnt1);
        end
    end
    agl = temp;
end

