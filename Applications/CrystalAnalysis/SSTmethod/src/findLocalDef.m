function [G,curlG] = findLocalDef(fff,angles,radii,masses,weights,numIter1,numIter2,plotFIG,GPUflag)
%% [G,curlG] = findLocalDef(angles,radii,masses,weights,numIter1,numIter2)
% given the result from the synchrosqueezed wave packet transform applied to an image of a hexagonal crystal,
% the code variationally identifies the local crystal deformation and outputs it as a matrix G for each pixel,
% where G is the deformation gradient of a deformation which maps the crystal onto a perfect, undistorted reference crystal
%
% the result from the synchrosqueezed wave packet transform is a measure;
% at each pixel this measure is the sum of three delta peaks in Fourier space;
% the position of those three peaks is described in polar coordinates via "angles" and "radii" and the mass of each peak via "masses" (each being m by n by 3)
%
% input:
% fff     - image
% angles  - polar angles of the synchrosqueezed wave packet transform peaks
% radii   - polar radii of the synchrosqueezed wave packet transform peaks
% masses  - masses of the synchrosqueezed wave packet transform peaks
% weights - weight of fidelity, curl, H^1, and elastic term (assuming grid spacing of 1)
% numIter1- number of Bregman iterations
% numIter2- number of additional Bregman iterations for the point group correction
% plotFIG - true if outputing figures
% GPUflag - true if computation on a GPU
%
% output:
% G       - the local deformation gradient at each pixel (m by n by 2 by 2)
% curlG   - the curl of the deformation gradient modulo the crystal point group

%% check input and set default parameters and flags
if nargin < 9
    GPUflag = false;
end
if nargin < 8
    plotFIG = false;
end
if nargin <7
    numIter2 = 20;
end
if nargin < 6
    numIter1 = 20;
end
if nargin < 5
    weights = [1 2 .1 .1];
end
lambda = 500;             % Bregman penalty weight
isLeastSquare = true;     % true if use a least square method to determine G0; false if use an average method to determine G0

waveVecs = cat(4,radii.*cos(angles),radii.*sin(angles));

%% create a stencil (with three vectors since we only consider one half-plane) for the dual lattice
stencil = zeros(3,2); % stencil in Fourier domain
M = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
stencil(1,1) = mean(radii(:));
for i = 2:3
    stencil(i,:) = stencil(i-1,:) * M;
end

%% produce an initial guess for G (note: we need G^T stencil = waveVec)
if ~isLeastSquare
    angle = mean(angles,3) - pi/3;
    radius = mean(radii,3) / norm(stencil(1,:),2);
    G0 = zeros([size(angle) 2 2]);
    G0(:,:,1,1) =  radius.*cos(angle);
    G0(:,:,2,2) =  radius.*cos(angle);
    G0(:,:,1,2) =  radius.*sin(angle);
    G0(:,:,2,1) = -radius.*sin(angle);
else
    G0 = zeros([size(radii(:,:,1)) 2 2]);
    A = [stencil(1,1) 0 stencil(1,2) 0;0 stencil(1,1) 0 stencil(1,2); ...
        stencil(2,1) 0 stencil(2,2) 0;0 stencil(2,1) 0 stencil(2,2); ...
        stencil(3,1) 0 stencil(3,2) 0;0 stencil(3,1) 0 stencil(3,2)];
    for k = 1:size(G0,1)
        for l = 1:size(G0,2)
            aux = A\[radii(k,l,1)*[cos(angles(k,l,1));sin(angles(k,l,1))]; ...
                radii(k,l,2)*[cos(angles(k,l,2));sin(angles(k,l,2))]; ...
                radii(k,l,3)*[cos(angles(k,l,3));sin(angles(k,l,3))]];
            G0(k,l,:) = aux([1 3 2 4]);
        end
    end
    G0(:,:,1,2) = -G0(:,:,1,2); % flipping direction of y-coordinate
    G0(:,:,2,1) = -G0(:,:,2,1);
end

%% visualize initial guess
if plotFIG
    [G0,curlG0] = BregmanIter_FitL1Curl(waveVecs,masses,stencil,G0,weights,lambda,0,GPUflag);
    visualize(fff,G0,curlG0,masses,0.15);
end

%% identify the local crystal deformation
if GPUflag
    gpuDevice(1);
end
tic;
[G1,curlG1] = BregmanIter_FitL1Curl(waveVecs,masses,stencil,G0,weights,lambda,numIter1,GPUflag);
toc

%% identify the local crystal deformation with wave vectors in [pi/6,7pi/6]
angles  = mod(angles-pi/6,pi)+pi/6;
angles2 = circshift(angles,[0 0 -1]);
radii2   = circshift(radii,[0 0 -1]);
aux  = repmat(angles2(:,:,3) > angles(:,:,3),[1 1 3]);
angles(aux) = angles2(aux);
radii(aux) = radii2(aux);
waveVecs2 = cat(4,radii.*cos(angles),radii.*sin(angles));
if ~isLeastSquare
    angle = mean(angles,3) - pi/3;
    G0(:,:,1,1) =  radius.*cos(angle);
    G0(:,:,2,2) =  radius.*cos(angle);
    G0(:,:,1,2) =  radius.*sin(angle);
    G0(:,:,2,1) = -radius.*sin(angle);
else
    G0 = zeros([size(radii(:,:,1)) 2 2]);
    for k = 1:size(G0,1)
        for l = 1:size(G0,2)
            aux = A\[radii(k,l,1)*[cos(angles(k,l,1));sin(angles(k,l,1))]; ...
                radii(k,l,2)*[cos(angles(k,l,2));sin(angles(k,l,2))]; ...
                radii(k,l,3)*[cos(angles(k,l,3));sin(angles(k,l,3))]];
            G0(k,l,:) = aux([1 3 2 4]);
        end
    end
    G0(:,:,1,2) = -G0(:,:,1,2); % flipping direction of y-coordinate
    G0(:,:,2,1) = -G0(:,:,2,1);
end
tic;
[G2,curlG2] = BregmanIter_FitL1Curl(waveVecs2,masses,stencil,G0,weights,lambda,numIter1,GPUflag);
toc

%% identify region where first segmentation has high curl, but second does not
gaussFilter = fspecial( 'gaussian', [5 5], 100 );
curl1 = imfilter(sqrt(sum(curlG1.^2,3)), gaussFilter, 'replicate');
curl2 = imfilter(sqrt(sum(curlG2.^2,3)), gaussFilter, 'replicate');
diffRegion = curl1-curl2 > .01;
diffRegion = imdilate(diffRegion,strel('disk',5));
diffRegion(curl2-curl1 > .01) = false;

%% perform segmentation with new rotation values
% Use G1 where first segmentation has low curl, G2 elsewhere
G3 = G1;
G3(repmat(diffRegion,[1 1 2 2])) = G2(repmat(diffRegion,[1 1 2 2]));
waveVecs(repmat(diffRegion,[1 1 3 2])) = waveVecs2(repmat(diffRegion,[1 1 3 2]));
% set rotation values at the boundaries of diffRegion
R = getRFromFinv(G3,6);
if 1
    tic;
    [G,curlG] = BregmanIter_FitL1Curl_PointGroup(waveVecs,masses,stencil,R,G3,weights,lambda,numIter2,GPUflag);
    toc
end

%% visualize the result
if plotFIG
    visualizeSST(fff,G,curlG,masses,0.05);
end
