function G0 = initGFromSSPeaks( angles, radii, stencilAngles, stencilRadii )
% given the positions of the Fourier peaks of the perfect lattice (stencilAngles,stencilRadii),
% compute deformation gradient matrices G0 at each pixel transforming the image lattice into the perfect lattice
% (or equivalently, G0^T transforms stencil into (angles,radii))
%
% angles,radii assume x- and y-coordinate to point into direction of increasing second and first index, respectively
% G0 assumes x- and y-coordinate to point into direction of increasing second and decreasing first index, respectively

% compute Fourier stencil of perfect lattice
stencil = [stencilRadii.*cos(stencilAngles) stencilRadii.*sin(stencilAngles)];

% set up overdetermined system to find G0
A = zeros(2*size(stencil,1),2);
A(1:2:end,1) = stencil(:,1);
A(1:2:end,3) = stencil(:,2);
A(2:2:end,2) = stencil(:,1);
A(2:2:end,4) = stencil(:,2);

% find G0
G0 = zeros([size(radii(:,:,1)) 2 2]);
for k = 1:size(G0,1)
    for l = 1:size(G0,2)
        rhs = permute([radii(k,l,:).*cos(angles(k,l,:)) radii(k,l,:).*sin(angles(k,l,:))],[2 1 3]); rhs = rhs(:);
        aux = A \ rhs;
        G0(k,l,:) = aux([1 3 2 4]);
    end
end

% flip direction of y-coordinate
G0(:,:,1,2) = -G0(:,:,1,2);
G0(:,:,2,1) = -G0(:,:,2,1);
end