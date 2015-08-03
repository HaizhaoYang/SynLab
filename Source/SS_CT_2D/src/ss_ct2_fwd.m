function [ss_energy coefCell coefTensor LocWavVecx LocWavVecy] = ss_ct2_fwd(img,SPg,R_low,R_high,NB,rad,is_real,is_unif,typeNUFFT,xo,yo,epsl,t_sc,s_sc,red,is_cos,wedge_length_coarse)
% ss_ct2_fwd.m - 2D Synchrosqueezed Curvelet Transform
%
% Inputs:
%   img         given image of size Nx by Ny
%
% Optional inputs:
%   SPg         SPg(1) is the number of samples in the vertical direction in the image
%               SPg(2) is the number of samples in the horizontal direction in the image
%   R_low       lower bound of interested frequency, not necessary in this code
%   R_high      upper bound of interested frequency
%   NB          (2*NB+1) by (2*NB+1) is the size of samples in the phase in [-R_high,R_high]*[-R_high,R_high]
%   rad         a parameter to adjust the size of supports of curvelets in
%               the frequency domain, rad <= 2.
%   is_real     0: img is complex
%               1: img is real
%   is_unif
%   typeNUFFT   1: NUFFT by Air Force Lab
%               2: USFFT by E. Candes
%               3: NUFFT by L. Greengard and J.-Y. Lee
%               4: Direct non-uniform Fourier Transform
%   xo          a matrix of size Nx by Ny, non-uniform locations in x at which img is measured
%   yo          a matrix of size Nx by Ny, non-uniform locations in y at which img is measured
%   epsl        threshold for synchrosqueezed curvelet transform
%               Note: The accuracy of the synchrosqueezed transform is determined by the
%               parameter epsl. However, due to numerical issues, epsl
%               cannot be too small.
%   t_sc    scaling parameter for radius
%               [default set to 1-1/8]
%   s_sc    scaling parameter for angle
%               [default set to 1/2+1/8]
%   red         a parameter for redundancy
%               [ default set to 1]
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   wedge_length_coarse
%               length of coarsest wedge
%               [default set to 8]
%
% Outputs:
%   ss_energy          synchrosqueezed energy distribution
%   coefCell           a cell structure storing general curvelet transform coefficients
%   coefTensor         a tensor storing the general curvelet transform coefficients
%   LocWavVecx         local wave vector in the first direction
%   LocWavVecy         local wave vector in the second direction
%
%by Haizhao Yang

% TODO: insert non-uniform transform

[Nx Ny] = size(img);
N = max(Nx,Ny);
if nargin < 2, SPg = round([Nx Ny]/8); end;
%if nargin < 3, R_low = 0; end;
if nargin < 4, R_high = ceil(N/2); end;
if nargin < 5, NB = ceil(R_high/8); end;
if nargin < 6, rad = 1; end;
if nargin < 7, is_real = 0; end;
if nargin < 8, is_unif = 1; end;
if nargin < 9, typeNUFFT = 1; end;
if nargin < 10, xo = 0:1/Nx:(Nx-1)/Nx; yo = 0:1/Ny:(Ny-1)/Ny; [xo yo] = ndgrid(xo,yo); end;
if nargin < 11, yo = repmat(0:1/Ny:(Ny-1)/Ny,Nx,1); end
if nargin < 12, epsl = 1e-2; end;
if nargin < 13, t_sc = 1 - 1/8; end;
if nargin < 14, s_sc = 1/2 + 1/8; end;
if nargin < 15, red = 1; end;
if nargin < 16, is_cos = 1; end;
if nargin < 17, wedge_length_coarse = 8; end;

R_low = 0;
[coefCell aaa bbb] = gdct2_fwd_red(img,is_real,SPg,R_high*sqrt(2),R_low,rad,is_cos, t_sc, s_sc, red,wedge_length_coarse);

ncl = numel(coefCell);
[t1 t2] = size(coefCell{1}{1});
nclp = 0;
for g = 1:ncl
    nclp = nclp + numel(coefCell{g});
end
coefTensor = zeros(t1,t2,nclp);
aa = zeros(t1,t2,nclp);
bb = zeros(t1,t2,nclp);
cnt_nclp = 1;
for g=1:ncl
    [szmccc,sznccc] = size(coefCell{g});
    for cnt1 = 1:szmccc
        for cnt2 = 1:sznccc
            coefTensor(:,:,cnt_nclp) = coefCell{g}{cnt1,cnt2};
            aa(:,:,cnt_nclp) = aaa{g}{cnt1,cnt2};
            bb(:,:,cnt_nclp) = bbb{g}{cnt1,cnt2};
            cnt_nclp = cnt_nclp + 1;
        end
    end
end
EXT = 10^10;
szc = size(coefTensor);
gud = find(abs(coefTensor)>epsl);
aa = real( aa(gud)./coefTensor(gud) / (2*pi*i));
bb = real( bb(gud)./coefTensor(gud) / (2*pi*i));
%dist = sqrt(aa.^2+bb.^2);
%good = find(dist<=R_high & dist>R_low & abs(aa)<Nx/2 & abs(bb)<Ny/2);
good = find(abs(aa)<=min(R_high,Nx/2) & abs(bb)<=min(R_high,Ny/2));

LocWavVecx = repmat(EXT,szc);
LocWavVecx(gud(good)) = aa(good);
LocWavVecy = LocWavVecx;
LocWavVecy(gud(good)) = bb(good);

if(1)
    WB = NB/R_high;
    ss_energy = zeros(2*NB+1,2*NB+1,t1,t2);
    for a=1:t1
        for b=1:t2
            tc = coefTensor(a,b,:);
            tk1 = LocWavVecx(a,b,:);
            tk2 = LocWavVecy(a,b,:);
            tc = tc(:);
            gud = find(tk1<EXT);
            tk1 = tk1(gud);
            tk2 = tk2(gud);
            for g=1:length(gud)
                loc1 = round(tk1(g)*WB)+NB+1;
                loc2 = round(tk2(g)*WB)+NB+1;
                ss_energy(loc1,loc2,a,b) = ss_energy(loc1,loc2,a,b) + abs(tc(gud(g))).^2;
            end
        end
    end
end


