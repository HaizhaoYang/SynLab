function [aaa,bbb,ccc] = wp3_ext(Nx,Ny,Nz,SPx,SPy,SPz,dat,rad, xxx,is_real,R_low,R_high)
% This code generate the partial derivative of the coefficients of the T3D
% wave packet transform with a fixed scaling parameter 1/2.
%
% Input:
% [Nx Ny Nz] is the size of the image
% [SPx SPy SPz] is the size of the samples in space
% dat stores the information of the grid of wave packets
% rad is the smallest size of the supports of wave packets
% xxx: given signal
% is_real: if xxx is real, is_real = 1, otherwise is_real = 0
% [R_low R_high]         The range of interested spectrum
%
% Output:
% aaa, bbb, ccc, derivatives in three directions
%
% By Haizhao Yang

N = max([Nx,Ny,Nz]);
pat = freq_pat(N/2,rad,R_low,R_high);

fff = fftn(xxx)/sqrt(numel(xxx));

%enlarge fff
L = floor(Nx/2);
fff = [fff(1:L,:,:); zeros(N-Nx,Ny,Nz); fff(L+1:end,:,:)];
L = floor(Ny/2);
fff = [fff(:,1:L,:) zeros(N,N-Ny,Nz) fff(:,L+1:end,:)];
temp = zeros(N,N,N);
L = floor(Nz/2);
temp(:,:,1:L) = fff(:,:,1:L);
temp(:,:,N-Nz+L+1:end) = fff(:,:,L+1:end);
fff = temp;

aaa = cell(numel(dat),1);
bbb = cell(numel(dat),1);
ccc = cell(numel(dat),1);

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
    
    t0 = fff(tmp1,tmp2,tmp3).*conj(cur.ex1); %local fourier
    t1 = zeros(SPx,SPy,SPz);
    t1(new1,new2,new3) = t0;
    aaa{I} = ifftn(t1)*sqrt(numel(t1));
    
    t0 = fff(tmp1,tmp2,tmp3).*conj(cur.ex2); %local fourier
    t1 = zeros(SPx,SPy,SPz);
    t1(new1,new2,new3) = t0;
    bbb{I} = ifftn(t1)*sqrt(numel(t1));
    
    t0 = fff(tmp1,tmp2,tmp3).*conj(cur.ex3); %local fourier
    t1 = zeros(SPx,SPy,SPz);
    t1(new1,new2,new3) = t0;
    ccc{I} = ifftn(t1)*sqrt(numel(t1));
end
end
