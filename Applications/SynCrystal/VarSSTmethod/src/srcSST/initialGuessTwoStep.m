function [radii angles BD TTEng_1st TTEng_2nd W_sec] = initialGuessTwoStep(fff,subsampleRate,plotFIG,parameters)
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
    parameters = struct('extendRatio',0.1,'threRatio',0.1,'epsFreqBand',1e-2,'numWave',3,'numGridPerSector',30,'red',[2 1],'rad',2,'s_sc',0.8,'t_sc',0.8,'isPeriodic',1,'epsSST',1e-2);
end
if nargin<3, plotFIG = false; end;
if nargin<2, subsampleRate = 2; end;
%% normalize the image
fff = fff/max(max(abs(fff)));

isPeriodic = parameters.isPeriodic;

%% estimate Frequency band
tic;
extendRatio = parameters.extendRatio;%influence the estimate around the boundaries.
threRatio = parameters.threRatio;
epsl = parameters.epsFreqBand;

if ~isPeriodic
    [m n] = size(fff);
    mx = m;
    nx = n;
    dm = round(m*0.2/subsampleRate)*subsampleRate;
    dn = round(n*0.2/subsampleRate)*subsampleRate;
    fff = [zeros(dm,dn*2+n);zeros(m,dn) fff zeros(m,dn);zeros(dm,dn*2+n)];
else
    [mx nx] = size(fff);
    dm = 0; dn = 0;
end

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
%set up parameters for SST
SPg = round(size(X)/subsampleRate);
R_low = st; R_high = ed;
num_wave = parameters.numWave;
NB = [(R_high - R_low),num_wave*parameters.numGridPerSector];
red = parameters.red;
rad = parameters.rad;
epsl = parameters.epsSST;
is_real = parameters.is_real;
t_sc = parameters.t_sc;
s_sc = parameters.s_sc;
num_direction = 1;
%----------------------------------------------------------------------
%ss transform
[ss_energy ss_avgdx ss_avgdy] = SS_ct_polar_v2(num_direction,X,SPg,NB,rad,is_real,[R_low,R_high],[0 2*pi],epsl,red,t_sc, s_sc);
ss_energy = permute(ss_energy,[3 4 1 2]);
if type == 0  
    ss_energy = ss_energy(:,1:round(nnx/subsampleRate),:,:);
else
    ss_energy = ss_energy(1:round(mmx/subsampleRate),:,:,:);
end
ss_energy = ss_energy(dm/subsampleRate+1:dm/subsampleRate+round(mx/subsampleRate),dn/subsampleRate+1:dn/subsampleRate+round(nx/subsampleRate),:,:);
ss_energy = permute(ss_energy,[3 4 1 2]);
ss_avgdx = permute(ss_avgdx,[3 4 1 2]);
if type == 0  
    ss_avgdx = ss_avgdx(:,1:round(nnx/subsampleRate),:,:);
else
    ss_avgdx = ss_avgdx(1:round(mmx/subsampleRate),:,:,:);
end
ss_avgdx = ss_avgdx(dm/subsampleRate+1:dm/subsampleRate+round(mx/subsampleRate),dn/subsampleRate+1:dn/subsampleRate+round(nx/subsampleRate),:,:);
ss_avgdx = permute(ss_avgdx,[3 4 1 2]);
ss_avgdy = permute(ss_avgdy,[3 4 1 2]);
if type == 0  
    ss_avgdy = ss_avgdy(:,1:round(nnx/subsampleRate),:,:);
else
    ss_avgdy = ss_avgdy(1:round(mmx/subsampleRate),:,:,:);
end
ss_avgdy = ss_avgdy(dm/subsampleRate+1:dm/subsampleRate+round(mx/subsampleRate),dn/subsampleRate+1:dn/subsampleRate+round(nx/subsampleRate),:,:);
ss_avgdy = permute(ss_avgdy,[3 4 1 2]);

%first step angle estimate
[aglPos TTEng_1st TTEng_2nd] = LocSmooth(sum(ss_energy),num_wave);
BD = 1./sqrt(TTEng_1st-TTEng_2nd+1);

%----------------------------------------------------------------------
%adjust the size of intAgl
[~,~,szm1,szn1] = size(ss_energy);
[szm2 szn2] = size(aglPos);
[gdx1 gdy1] = meshgrid(0:1/(szm1-1):1,0:1/(szn1-1):1);
[gdx2 gdy2] = meshgrid(0:1/(szm2-1):1,0:1/(szn2-1):1);
aglPos = (interp2(gdx2,gdy2,aglPos',gdx1,gdy1,'nearest'))';
BD = (interp2(gdx2,gdy2,BD',gdx1,gdy1))'; 
intAgl= aglPos*180/NB(2);
if plotFIG
    figure;imagesc(intAgl);axis image;colorbar;
end
aglPos = round(aglPos);
%----------------------------------------------------------------------
%second step angle estimate and bump detection
[angles radii TTEng_1st TTEng_2nd W_sec] = LocWeight_v2(ss_energy,ss_avgdx,ss_avgdy,num_wave,aglPos);

if parameters.isAffine == 1
    radii = radii/(size(X,1));
end

%check aglPos
if 0
    NC = 20;
    LL = round(parameters.numGridPerSector/4);
    [m n] = size(aglPos);
    pos1 = round(rand(1,NC)*(m-NC))+(1:NC);
    pos2 = round(rand(1,NC)*(n-NC))+(1:NC);
    figure;
    for cnt = 1:NC
        temp = ss_energy(:,:,pos1(cnt),pos2(cnt));
        subplot(4,1,1);imagesc(temp);axis image;
        temp2 = zeros(size(temp));
        for cnt2 = 1:num_wave
            cnt3 = -LL:LL;
            d = cnt3 + aglPos(pos1(cnt),pos2(cnt)) + cnt2*(NB(2))/num_wave;
            d = mod(d-1,NB(2)+1)+1;
            temp2(:,d) = temp(:,d);
        end
        subplot(4,1,2);imagesc(temp2);axis image;
        subplot(4,1,3);imagesc(ones(NB(1),1)*checktemp(:,pos1(cnt),pos2(cnt))');axis image;
        temp2 = zeros(size(temp));
        d = aglPos(pos1(cnt),pos2(cnt)) + 0:(NB(2))/num_wave:NB(2);
        temp2(:,d) = 1;
        subplot(4,1,4);imagesc(temp2);axis image;
        head = sprintf('%d',aglPos(pos1(cnt),pos2(cnt)));
        title(head);
        pause;
    end
end

clear ss_energy ss_avgdx ss_avgdy;
%----------------------------------------------------------------------

