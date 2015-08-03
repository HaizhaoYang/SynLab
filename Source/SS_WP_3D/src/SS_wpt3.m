function [kb sss ss kk1 kk2 kk3] = SS_wpt3(fff,SPg,NB,rad,is_real,R_low,R_high,epsl)
% This code computes the 3D synchrosqueezed wave packet transform with a fixed
% scaling parameter 1/2.
%
% Input:
% fff is the image
% SPg(1) is the number of samples in the x direction in the image
% SPg(2) is the number of samples in the y direction in the image
% SPg(3) is the number of samples in the z direction in the image
% 2*NB(1)+1 by 2*NB(2)+1 bt 2*NB(3)+1 is the size of samples in the phase in
% [-R_high,R_high]*[-R_high,R_high]*[-R_high,R_high]
% rad is the smallest size of supports of wave packets
% is_real = 0: fff is complex
% is_real = 1: fff is real
% R_low: lower bound of interested frequency
% R_high: upper bound of interested frequency
% epsl: the threshold for small coefficients of wave packet transform
% Note: The accuracy of the synchrosqueezed transform is determined by the
% parameter epsl.
%
% Output:
% kb the synchrosqueezed energy distribution
% sss the wave packet transform coefficients
% ss a tensor storing the wave packet transform coefficients
%
% kk1: local wave vector in the first direction
% kk2: local wave vector in the second direction
% kk3: local wave vector in the third direction
%
% By Haizhao Yang

%set up size of image
[Nx Ny Nz] = size(fff);
N = max([Nx,Ny,Nz]);
%set up size of samples in the space
SPx = SPg(1);
SPy = SPg(2);
SPz = SPg(3);

dat = wp3_setup(Nx,Ny,Nz,SPx,SPy,SPz,rad,R_low,R_high);
%sss is the wave atom coefficients
sss = wp3_fwd(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, fff, is_real,R_low,R_high);
%aaa and bbb are the wave atom coefficients with the derivatives in b_1 and b_2
[aaa,bbb,ccc] = wp3_ext(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, fff, is_real,R_low,R_high);
ncl = numel(sss);
ss = zeros(SPx,SPy,SPz,ncl);
aa = zeros(SPx,SPy,SPz,ncl);
bb = zeros(SPx,SPy,SPz,ncl);
cc = zeros(SPx,SPy,SPz,ncl);
for g=1:ncl
    ss(:,:,:,g) = sss{g};
    aa(:,:,:,g) = aaa{g};
    bb(:,:,:,g) = bbb{g};
    cc(:,:,:,g) = ccc{g};
end
clear aaa bbb ccc;
gud = find(abs(ss)>epsl);
aa = round(real( aa(gud)./ss(gud) / (2*pi*i)));
bb = round(real( bb(gud)./ss(gud) / (2*pi*i)));
cc = round(real( cc(gud)./ss(gud) / (2*pi*i)));
EXT = 10^10;
szs = size(ss);
dist = sqrt(aa.^2+bb.^2+cc.^2);
good = find(dist<=R_high & dist>R_low & abs(aa)<Nx/2 & abs(bb)<Ny/2 & abs(cc)<Nz/2);
kk1 = repmat(EXT,szs);
kk1(gud(good)) = aa(good);
kk2 = kk1;
kk2(gud(good)) = bb(good);
kk3 = kk1;
kk3(gud(good)) = cc(good);
clear dist;

WB = NB/R_high;
kb = zeros(2*NB(1)+1,2*NB(1)+1,2*NB(3)+1,SPx,SPy,SPz);
for a=1:SPx
    for b=1:SPy
        for c=1:SPz
            tc = ss(a,b,c,:);
            tk1 = kk1(a,b,c,:);
            tk2 = kk2(a,b,c,:);
            tk3 = kk3(a,b,c,:);
            tc = tc(:);
            gud = find(tk1<EXT);
            tk1 = tk1(gud);
            tk2 = tk2(gud);
            tk3 = tk3(gud);
            for g=1:length(gud)
                loc1 = round( tk1(g)*WB(1)) + NB(1) + 1;
                loc2 = round( tk2(g)*WB(2)) + NB(2) + 1;
                loc3 = round( tk3(g)*WB(3)) + NB(3) + 1;
                kb(loc1,loc2,loc3,a,b,c) = kb(loc1,loc2,loc3,a,b,c) + abs(tc(gud(g))).^2;
            end
        end
    end
end


