function fff = SS_wpt3_inv(kb,sss,ss,kk1,kk2,kk3,sz,NB,rad,R_low,R_high,is_real)
% This code computes the 3D inverse synchrosqueezed wave packet transform
% with a fixed scaling parameter 1/2
%
% Input:
% kb: the synchrosqueezed energy distribution
% sss: the wave packet transform coefficients stored in cells
% ss: the wave packet transform coefficients stored in a tensor
% kk1: local wave vector in the x direction
% kk2: local wave vector in the y direction
% kk3: local wave vector in the z direction
% sz is the size of the original image
% 2*NB(1)+1 by 2*NB(2)+1 by 2*NB(3)+1: the size of grid in the phase domain
% in [-R_high,R_high]*[-R_high,R_high]*[-R_high,R_high]
% rad is the smallest size of supports of wave packets
% R_low: lower bound of interested frequency
% R_high: upper bound of interested frequency
% is_real = 0: fff is complex
% is_real = 1: fff is real
%
% Output:
% fff: the original image
%
% By Haizhao Yang

Nx = sz(1);
Ny = sz(2);
Nz = sz(3);
N = max([Nx,Ny,Nz]);
dd = zeros(size(ss));
[SPx SPy SPz] = size(sss{1});
dat = wp3_setup(Nx,Ny,Nz,SPx,SPy,SPz,rad,R_low,R_high);
EXT = 10^10;
ncl = numel(sss);
WB = NB/R_high;
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
                if( kb(loc1,loc2,loc3,a,b,c)>0)
                     dd(a,b,c,gud(g)) = ss(a,b,c,gud(g));
                end
            end
        end
    end
                    
end
ddd = cell(size(sss));
for g=1:ncl
    ddd{g} = dd(:,:,:,g);
end
fff = wp3_inv(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, ddd,is_real,R_low,R_high);