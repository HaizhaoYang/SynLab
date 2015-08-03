function [kb kk1 kk2] = SS_ct_polar(num_direction,fff,SPg,NB,rad,is_real,R_low,R_high,epsl,red,t_sc, s_sc,is_cos,  wedge_length_coarse)
% Highly redundant SST
% Generate multiframes by rotation in the frequency domain
%
% Input:
% fff is the image
% SPg(1) is the number of samples in the vertical direction in the image
% SPg(2) is the number of samples in the horizontal direction in the image
% NB(1) is the size of radial grid in [R_low R_high]
% NB(2) is the size of angular grid in [0,pi)
% rad is the smallest size of supports of wave packets
% is_real = 0: fff is complex
% is_real = 1: fff is real
% R_low: lower bound of interested frequency
% R_high: upper bound of interested frequency
% epsl: the threshold for small coefficients of general curvelet transform
% Note: The accuracy of the synchrosqueezed transform is determined by the
% parameter epsl.
% red is the redundancy parameter, default set to 1
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   t_sc    scaling parameter for radius
%               [default set to 1-1/8]
%   s_sc    scaling parameter for angle
%               [default set to 1/2+1/8]
%   wedge_length_coarse
%               length of coarsest wedge
%               [default set to 5]
%
% Output:
% kb the synchrosqueezed energy distribution in a polar coordinate
%   size(kb) = NB,
% ccc the general curvelet transform coefficients
% cc a tensor storing the general curvelet transform coefficients
%
% kk1: local wave vector in the first direction
% kk2: local wave vector in the second direction
%
%by Haizhao Yang

if nargin < 10, red = 1; end;
if nargin < 11, t_sc = 1 - 1/8; end;
if nargin < 12, s_sc = 1/2 + 1/8; end;
if nargin < 13, is_cos = 1; end;
if nargin < 14, wedge_length_coarse = 8; end;

%set up size of image
[Nx Ny] = size(fff);
N = max(Nx,Ny);
%set up size of samples in the space
SPx = SPg(1); SPy = SPg(2);

%ccc is the general curvelet coefficients
%aaa and bbb are the general curvelet coefficients with the derivatives in b_1 and b_2
[ccc aaa bbb] = gdct2_fwd_red(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos, t_sc, s_sc, red, wedge_length_coarse);

ncl = numel(ccc);
[t1 t2] = size(ccc{1}{1,1});
nclp = 0;
for g = 1:ncl
    nclp = nclp + numel(ccc{g});
end
cc = zeros(t1,t2,nclp);
aa = zeros(t1,t2,nclp);
bb = zeros(t1,t2,nclp);
cnt_nclp = 1;
for g=1:ncl
    [szmccc,sznccc] = size(ccc{g});
    for cnt1 = 1:szmccc
        for cnt2 = 1:sznccc
            cc(:,:,cnt_nclp) = ccc{g}{cnt1,cnt2};
            aa(:,:,cnt_nclp) = aaa{g}{cnt1,cnt2};
            bb(:,:,cnt_nclp) = bbb{g}{cnt1,cnt2};
            cnt_nclp = cnt_nclp + 1;
        end
    end
end
clear aaa bbb ccc;
EXT = 10^10;
szc = size(cc);
% kk1 local wave vector in the first direction
% kk2 local wave vector in the second direction
gud = find(abs(cc)>epsl);
aa = round(real( aa(gud)./cc(gud) / (2*pi*i)));
bb = round(real( bb(gud)./cc(gud) / (2*pi*i)));
dist = sqrt(aa.^2+bb.^2);
good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);

kk1 = repmat(EXT,szc);
kk1(gud(good)) = aa(good);
kk2 = kk1;
kk2(gud(good)) = bb(good);
clear aa bb;

%temp = kk1(gud(good));
%median(temp(:))
%temp = kk2(gud(good));
%median(temp(:))

if(1)
    dr = (R_high-R_low)/(NB(1));
    da = pi/(NB(2))/num_direction;
    %WB = NB/R_high;%2*min(N,2*R_high)/NB;
    kb = zeros(NB(1)+1,NB(2)+1,t1,t2);
    kb = SS_polar_old(cc,kk1,kk2,EXT,pi/num_direction,da,dr,NB(1)+1,NB(2)+1,R_low);
end


