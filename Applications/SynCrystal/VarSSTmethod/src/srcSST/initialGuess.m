function [R agl TTEng_1st TTEng_2nd W_sec ss_energy] = initialGuess(fff,subsampleRate,plotFIG,parameters)
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

if 0
    %rotate 90 degrees
    [m n] = size(fff);
    temp = zeros(n,m);
    for cnt1 = 1:m
        for cnt2 = 1:n
            temp(cnt1,cnt2) = fff(n+1-cnt2,cnt1);
        end
    end
    fff = temp;
end
%% estimate Frequency band
tic;
extendRatio = parameters.extendRatio;%influence the estimate around the boundaries.
threRatio = parameters.threRatio;
epsl = parameters.epsFreqBand;
[m n] = size(fff);
if m ~= n
    mm = max(m,n);
    if m > n
        X = [fff zeros(m,m-n)];
    else
        X = [fff; zeros(n-m,n)];
    end
    idx = 1:mm;
    X = X(idx,idx);
else
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
%set up parameters for SST
SPg = round(size(fff)/subsampleRate);
R_low = st; R_high = ed;
%set up the size of samples in [-R_high,R_high]*[-R_high,R_high] in the
%phase domain
num_wave = parameters.numWave;%determined by the reference grain
NB = [(R_high - R_low),num_wave*parameters.numGridPerSector];
%set up parameters for ss transform
red = parameters.red;
rad = parameters.rad;
epsl = parameters.epsSST;
is_real = parameters.is_real;
t_sc = parameters.t_sc;
s_sc = parameters.s_sc;
num_direction = 1;
dm = 0; dn = 0;
%----------------------------------------------------------------------
%ss transform
[ss_energy ss_avgdx ss_avgdy] = SS_ct_polar(num_direction,fff,SPg,NB,rad,is_real,R_low,R_high,epsl,red,t_sc, s_sc);

%----------------------------------------------------------------------
%bump detection
[agl R TTEng_1st TTEng_2nd W_sec] = LocWeight(ss_energy,ss_avgdx,ss_avgdy,num_wave);
clear ss_avgdx ss_avgdy;
%----------------------------------------------------------------------

