% This code uses the second method in "Crystal image analysis via 2D
% synchrosqueezed tranforms" to extract grain boundary and orientation.
% Based on this information, a variational optimization is applied to obtain
% results with more physical meaning. The variational model is detailed in
% "Combining 2D synchrosqueezed wave packet transform with optimization for
% crystal image analysis".
%
% By Jianfeng Lu, Benedikt Wirth and Haizhao Yang

function [phi,G,curlG,G0,curlG0,binaryIndex]=varSSTns(example,IsFC)
%% -------------------------------------------------------------------------
% Inputs:
% example       integer from 1 to 6 corresponding to 6 different crystal
%               images
% IsFC          0: coarseThre = fineThre, using the same defect region for
%                  computing Burger's vector and optimization
%               1: coarseThre < fineThre, using coarse and fine defect region
%                  for computing Burger's vector and optimization,
%                  respectively
%
% Outputs:
% phi           original image
% G             optimized gradient field
% curlG         curl of G
% G0            initial gradient field provided by SST
% binaryIndex   fine defect regions with indices identifying connected
%               components
%% -------------------------------------------------------------------------
if nargin <1, example = 1; IsFC = 1; end;
if nargin <2, IsFC = 1; end;

%% -------------------------------------------------------------------------
% Generate data and choose parameters
%
% 1) parameters for frequency band identification
% 1.1) 'extendRatio', usually takes 0 to 0.1, will influence the radii
%      estimation on the boundary but not sensitive
% 1.2) 'threRatio', a parameter to combine several possible close energy
%      peaks in the frequency band identification; threRatio=0.1 by
%      default; results are not sensitive to it.
% 1.3) 'epsFreqBand', denoising threshold when identifying the frequency
%      band; epsFreqBand=1e-2 by default
%
% 2) parameters for the synchrosqueezed transform
% 2.1) 'is_real',0,'numWave', the number of local wave vectors of interest, determined
%      by the crystal, e.g., numWave=3 for hexagonal crystals
% 2.2) 'numGridPerSector', the number of grid points in one sector;
%      in the case of hexagonal crystals, each sector covers an angle of
%      pi/3; numGridPerSector = 30 by default; larger number will not give
%      more accurate results
% 2.3) 'red', a redundancy parameter of the wave packet transform; red=1 by
%      default; the larger the more robust to noise, but the code is slower
% 2.4) 'rad' control the size of the support of the mother wave packet in
%      the K-plane; rad in [0.5,2]; larger rad gives finer results while
%      smaller rad gives smoother and more robust results;
% 2.5) 's_sc' and 't_sc' are scaling parameters for wave packet transforms;
%       usually 0.8 <=s_sc <= t_sc <= 1 for clean data; 0.625<=s_sc <=0.75,
%       t_sc = 0.8 for noisy data;
% 2.6) 'epsSST', denoising threshold; epsSST=1e-2 by default
%
% 3) parameters for defect region identification
% 3.1) 'energyThre', a thresholding parameter to identify valleys from
%      energy landmass provided by SST; energyThre in [0,3], smaller vthre
%      leads to more narrow boundaries;
% 3.2) 'coarseThre', a thresholding parameter to generate a coarse binary
%      image identifying grains and grain boundaries;
%      coarseThre in [0,1]; smaller value leads to wider boundaries
% 3.3) 'fineThre', a thresholding parameter to generate a fine binary
%      image identifying grains and grain boundaries;
%      fineThre in [0,1]; larger value leads to more narrow boundaries
%
% 4) parameters for the optimization
% 4.1) 'numIter', number of interations;
% 4.2) 'disVolThre', a thresholding parameter to show the distortion
%      volumne;
%% -------------------------------------------------------------------------
is_Bin = 0;
switch example
    case 1 
        % A PFC image
        load 'GB1.mat';
        % set up noise level
        NM = 0.5;%if NM = 1, then SRN = -8.34; if NM = 0.5, then SNR = -2.35;
        if NM>1
            parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                'epsFreqBand',1e-2,'is_real',0,'numWave',3, ...
                'numGridPerSector',30,'red',8,'rad',2, ...
                's_sc',0.6,'t_sc',0.6,'epsSST',0.15,'is_Bin',is_Bin);
            subsampleRate = 4;
        else if NM>0.5
                parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                    'epsFreqBand',1e-2,'is_real',0,'numWave',3, ...
                    'numGridPerSector',30,'red',10,'rad',2, ...
                    's_sc',0.625,'t_sc',0.625,'epsSST',0.15,'is_Bin',is_Bin);
                subsampleRate = 4;% NM <= 1
            else
                parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                    'epsFreqBand',1e-2,'is_real',0,'numWave',3, ...
                    'numGridPerSector',30,'red',10,'rad',2, ...
                    's_sc',0.625,'t_sc',0.625,'epsSST',0.15,'is_Bin',is_Bin);
                subsampleRate = 4;% NM <= 1
            end
        end
        energyThre = 2.4;
        disVolThre = 0.05;
        numIter = 60;
        if ~IsFC
            coarseThre = 0.3; fineThre = 0.3;  
            fprintf('Using coarseThre = fineThre = %f\n',fineThre);
        else
            coarseThre = 0.2; fineThre = 0.4;
            fprintf('Using coarseThre %f;   fineThre = %f\n',coarseThre,fineThre);
        end
end

% normalize the given image so as to use uniform thresholding parameters
fff = phi;
fff = fff/max(abs(fff(:)));
fff = fff - mean(fff(:));
[Nx,Ny] = size(fff);
ns = real(NM*(randn(Nx,Ny)+i*randn(Nx,Ny)));
snr = 10*log10((var(fff(:)))/var(ns(:)))
fff = fff + ns;
figure;plot(phi(:,1));
figure;plot(fff(:,1));

%% -------------------------------------------------------------------------
% The synchrosqueezed wave packet transform provids an initial guess
tic;
plotFIG = false;
[radii angles TTEng_1st TTEng_2nd masses] = initialGuess(fff,subsampleRate,plotFIG,parameters); 
% assumes y-coordinate pointing into direction of growing first index

% bring wave vectors into order according to angle
[angles,perm] = sort(angles,3);
for k = 1:size(perm,1)
    for l = 1:size(perm,2)
        radii(k,l,:) = radii(k,l,perm(k,l,:));
    end
end

if 0
    % Show angles and radii
    figure;
    subplot(2,3,1);imagesc(angles(:,:,1));axis image;colorbar; title('new angle 1');
    subplot(2,3,2);imagesc(angles(:,:,2));axis image;colorbar; title('new angle 2');
    subplot(2,3,3);imagesc(angles(:,:,3));axis image;colorbar; title('new angle 3');
    subplot(2,3,4);imagesc(radii(:,:,1));axis image;colorbar; title('new radii 1');
    subplot(2,3,5);imagesc(radii(:,:,2));axis image;colorbar; title('new radii 2');
    subplot(2,3,6);imagesc(radii(:,:,3));axis image;colorbar; title('new radii 3');
end

%% -------------------------------------------------------------------------
% find defect region and their connected components
if 1
    binaryIndexCoarse = 1-findDeftArea(smoothImage( sum(masses,3), 4, 1),energyThre,coarseThre);
    [CCCoarse,colComponentsCoarse] = periodicConnComp( binaryIndexCoarse );
    binaryIndex = 1-findDeftArea(smoothImage( sum(masses,3), 4, 1),energyThre,fineThre);
end
[CC,colComponents] = periodicConnComp( binaryIndex );

if 0
    pic = figure; imagesc(colComponentsCoarse); axis image; title('coarse defect regions w/ p.b.c.'); colorbar;
    fileName = sprintf('../../results/Part/nsCoarseDefectEx%d',example);colormap (1-gray);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    saveas(pic,fileName);
    pic = figure; imagesc(colComponents); axis image; title('defect regions w/ p.b.c.'); colorbar;
    fileName = sprintf('../../results/Part/nsFineDefectEx%d',example);colormap (1-gray);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    saveas(pic,fileName);
end
toc

%% -------------------------------------------------------------------------
tic;
% find initial noisy deformation gradient, assuming y-coordinate pointing into direction of decreasing 1st index, and using two different point group interpretations
G1 = initGFromSSPeaks( angles, radii, pi/3*[0;1;2], median(radii(:))*ones(3,1) );

%% -------------------------------------------------------------------------
% find locations of point group effect (first define the point group)
rot = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
pointGroup = eye(2);
for j = 1:5
    pointGroup = cat(3,pointGroup,pointGroup(:,:,j)*rot);
end
[L,G0] = findPointGroupJumpsBySweeping( pointGroup, G1, colComponentsCoarse );

%% -------------------------------------------------------------------------
% compute curl G0 on each connected component
curlG0 = findCurl( G0, pointGroup, L );
stencil = zeros(3,2); % stencil in Fourier domain
M = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
stencil(1,1) = median(radii(:));
for ii = 2:3
    stencil(ii,:) = stencil(ii-1,:) * M;
end
latticeVectors = [cos(pi/6) -sin(pi/6);sin(pi/6) cos(pi/6)]*stencil(1:2,:)'/(stencil(1,1)^2+stencil(1,2)^2)*length(G0); % Bravais vectors, assuming pixel distance is 1 and that mother wave packet lies on k1-axis in Fourier domain
prescribedCurl = findBurgersVectors( curlG0, latticeVectors, CCCoarse, CC );
toc

%% -------------------------------------------------------------------------
% show initialization in the whole domain
address = sprintf('../../results/Part/ns');
fileName = sprintf('IntWhole');
visualize(phi,G0,curlG0,masses,disVolThre,binaryIndex,0,address,fileName,example, CC);

%% -------------------------------------------------------------------------
% elastic denoising
tic;
G = elasticDenoising( G0, binaryIndex, prescribedCurl, pointGroup, L, numIter );
curlG = findCurl( G, pointGroup, L );
toc

%% -------------------------------------------------------------------------
% show optimized results in the whole domain
address = sprintf('../../results/Part/ns');
fileName = sprintf('OptWhole');
visualize(fff,G,curlG,masses,disVolThre,binaryIndex,0,address,fileName,example, CC);


