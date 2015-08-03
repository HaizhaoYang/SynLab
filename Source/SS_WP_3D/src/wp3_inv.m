function yyy = wp3_inv(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, ccc,is_real,R_low,R_high)
% This code generates the 3D inverse wave packet transform with a fixed
% scaling parameter 1/2.
%
% Input:
% [Nx Ny Nz] is the size of the image
% [SPx SPy SPz] is the size of the samples in space
% dat stores the information of the grid of wave packets
% rad is the smallest size of the supports of wave packets
% ccc is the variable for wave packet coefficients
% is_real: if xxx is real, is_real = 1, otherwise is_real = 0
% [R_low R_high]         The range of interested spectrum
%
% Output:
% yyy, signal
%
% By Haizhao Yang

N = max([Nx,Ny,Nz]);
pat = freq_pat(N/2,rad,R_low,R_high);

fff = zeros(N,N,N);
for I=1:numel(dat)
    cur = dat{I};
    tmp1 = mod(cur.k1gd,N)+1;
    tmp2 = mod(cur.k2gd,N)+1;
    tmp3 = mod(cur.k3gd,N)+1;
    len1 = numel(cur.k1gd);
    len2 = numel(cur.k2gd);
    len3 = numel(cur.k3gd);
    new1 = mod(cur.k1gd,SPx)+1;
    new2 = mod(cur.k2gd,SPy)+1;
    new3 = mod(cur.k3gd,SPz)+1;
    
    t1 = fftn(ccc{I})/sqrt(numel(ccc{I}));
    
    t0 = t1(new1,new2,new3);
    fff(tmp1,tmp2,tmp3) = fff(tmp1,tmp2,tmp3) + t0.*cur.wgt;
end
%cut fff
L = floor(Nx/2);
fff = [fff(1:L,:,:); fff(N-Nx+L+1:end,:,:)];
L = floor(Ny/2);
fff = [fff(:,1:L,:) fff(:,N-Ny+L+1:end,:)];
L = floor(Nz/2);
temp = zeros(Nx,Ny,Nz);
temp(:,:,1:L) = fff(:,:,1:L);
temp(:,:,L+1:end) = fff(:,:,N-Nz+L+1:end);
fff = temp;


yyy = ifftn(fff)*sqrt(numel(fff));
if is_real, yyy = real(yyy); end;
end

