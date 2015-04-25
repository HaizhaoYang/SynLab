function R = getRFromFinv(G,nsymm)
% function R = getRGromFinv(G,nsymm)
% R is m x n x 2 (for G being m x n x 2 x 2), giving the best rotations
% to minimize ||G - R*G(x,y+1)||, ||G - R*G(x+1,y)|| in L2-norm,
% where R is one of nsymm-fold symmetries.

[m,n,p,q] = size(G);
R = zeros(m,n,2);

alpha = linspace(0,2*pi,nsymm+1);
alpha = alpha(1:end-1);

mismatch = zeros(m,n,length(alpha));

for ii = 1:length(alpha)
  mismatch(:,:,ii) = ((G(:,:,1,1)-( cos(alpha(ii))*G([2:end 1],:,1,1)-sin(alpha(ii))*G([2:end 1],:,2,1))).^2 + ...
                      (G(:,:,1,2)-( cos(alpha(ii))*G([2:end 1],:,1,2)-sin(alpha(ii))*G([2:end 1],:,2,2))).^2 + ...
                      (G(:,:,2,1)-( sin(alpha(ii))*G([2:end 1],:,1,1)+cos(alpha(ii))*G([2:end 1],:,2,1))).^2 + ...
                      (G(:,:,2,2)-( sin(alpha(ii))*G([2:end 1],:,1,2)+cos(alpha(ii))*G([2:end 1],:,2,2))).^2);
end
%mismatch(:,:,1) = mismatch(:,:,1) / 10;
[tr,I] = min(mismatch,[],3);
R(:,:,1) = alpha(I);

for ii = 1:length(alpha)
  mismatch(:,:,ii) = ((G(:,:,1,1)-( cos(alpha(ii))*G(:,[2:end 1],1,1)-sin(alpha(ii))*G(:,[2:end 1],2,1))).^2 + ...
                      (G(:,:,1,2)-( cos(alpha(ii))*G(:,[2:end 1],1,2)-sin(alpha(ii))*G(:,[2:end 1],2,2))).^2 + ...
                      (G(:,:,2,1)-( sin(alpha(ii))*G(:,[2:end 1],1,1)+cos(alpha(ii))*G(:,[2:end 1],2,1))).^2 + ...
                      (G(:,:,2,2)-( sin(alpha(ii))*G(:,[2:end 1],1,2)+cos(alpha(ii))*G(:,[2:end 1],2,2))).^2);
end
%mismatch(:,:,1) = mismatch(:,:,1) / 10;
[tr,I] = min(mismatch,[],3);
R(:,:,2) = alpha(I);
