function bw = findDeftArea(v,vthre,bwthre)
% Inputs:
% v          an image containing several 'tablelands' as grains and valleys
%            as grain boundaries
% vthre      a thresholding parameter to identify valleys from v
%            vthre in [0,3], smaller vthre leads to more narrow boundaries
%            by default vthre = 2.25
% bwthre     a thresholding parameter to generate a binary image
%            identifying grains and grain boundaries
%            bwthre in [0,1], smaller bwthre leads to wider boundaries
%            by default bwthre = 0.5
%
% Outputs: 
% bw         a binary image identifying grains and grain boundaries

if nargin < 3, bwthre = 0.5; end;
if nargin < 2, vthre = 2.25; end;

bw = zeros(size(v)); bw(v<vthre) = 1;
[m n] = size(bw);
bw = (3-sharpen(1-bw,20,3,1))/3;

[pos1 pos2] = find(bw>0);
pos = sub2ind([m,n],pos1,pos2);
tmp = bw(sub2ind([m,n],mod(pos1-2,m)+1,pos2)) + bw(sub2ind([m,n],mod(pos1,m)+1,pos2)) +...
    bw(sub2ind([m,n],pos1,mod(pos2-2,n)+1)) + bw(sub2ind([m,n],pos1,mod(pos2,n)+1)) +...
    bw(sub2ind([m,n],mod(pos1-2,m)+1,mod(pos2,n)+1)) + bw(sub2ind([m,n],mod(pos1,m)+1,mod(pos2,n)+1)) +...
    bw(sub2ind([m,n],mod(pos1-2,m)+1,mod(pos2-2,n)+1)) + bw(sub2ind([m,n],mod(pos1,m)+1,mod(pos2-2,n)+1));
bw(pos(tmp==0)) = 0;

pos = find(bw>0);
for cnt = 1:5
    bw(pos) = 1;
    pos2 = find(bw>bwthre);
    bw(pos2) = 1;
    bw = smoothImage(bw,3,1);
end
bw = im2bw(bw,bwthre/2);

if 1
    for cnt = 1:5
        
        [pos1 pos2] = find(bw>0);
        pos = sub2ind([m,n],pos1,pos2);
        tmp = bw(sub2ind([m,n],mod(pos1-2,m)+1,pos2)) + bw(sub2ind([m,n],mod(pos1,m)+1,pos2));
        bw(pos(tmp==0)) = 0;
        tmp = bw(sub2ind([m,n],pos1,mod(pos2-2,n)+1)) + bw(sub2ind([m,n],pos1,mod(pos2,n)+1));
        bw(pos(tmp==0)) = 0;
        
        [pos1 pos2] = find(bw==0);
        pos = sub2ind([m,n],pos1,pos2);
        tmp = bw(sub2ind([m,n],mod(pos1-2,m)+1,pos2)) + bw(sub2ind([m,n],mod(pos1,m)+1,pos2));
        bw(pos(tmp==2)) = 1;
        tmp = bw(sub2ind([m,n],pos1,mod(pos2-2,n)+1)) + bw(sub2ind([m,n],pos1,mod(pos2,n)+1));
        bw(pos(tmp==2)) = 1;
        
    end
end