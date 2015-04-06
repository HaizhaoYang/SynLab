% In the first step, this code uses the first method in "Crystal image analysis via 2D
% synchrosqueezed tranforms" to extract grain boundary and orientation.
% Based on this information, this code applies the second method in
% the same paper to estimate crystal deformation and finer results of grain
% orientation. Finally, a variational optimization is applied to obtain
% results with more physical meaning. The variational model is detailed in
% "Combining 2D synchrosqueezed wave packet transform with optimization for
% crystal image analysis".
%
% By Jianfeng Lu, Benedikt Wirth and Haizhao Yang


function [phi,G,curlG,G0,curlG0,binaryIndex]=varSSTTwoStep(example,IsFC)
%% -------------------------------------------------------------------------
% Inputs:
% example       integer from 1 to 14
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
if nargin <1, example = 1; IsFC = 0; end;
if nargin <2, IsFC = 0; end;

switch example
    case 1
        phi = double(imread('GB13.jpg'));
        phi = double(phi(1:424,1:632));
        parameters = struct('extendRatio',0.05,'threRatio',0.05, ...
            'epsFreqBand',1e-3,'BD_type',1,'isAffine',0,'is_real',0,'numWave',3, ...
            'numGridPerSector',30,'red',[10 15],'rad',1, ...
            's_sc',0.65,'t_sc',0.75,'isPeriodic',0,'epsSST',1e-3);
        parameters2 = struct('extendRatio',0.05,'threRatio',0.05, ...
            'epsFreqBand',1e-3,'BD_type',1,'isAffine',1,'is_real',1,'numWave',3, ...
            'numGridPerSector',7,'red',[10 8],'rad',1, ...
            's_sc',0.6,'t_sc',0.75,'isPeriodic',0,'epsSST',1e-3);
        subsampleRate = 4;
        energyThre = 1.2;
        disVolThre = 0.1;
        BDThre = 0.95;
        numIter = 20;
        if ~IsFC
            coarseThre = 0.3; fineThre = 0.3;
            fprintf('Using coarseThre = fineThre = %f\n',fineThre);
        else
            coarseThre = 0.3; fineThre = 0.5;
            fprintf('Using coarseThre %f;   fineThre = %f\n',coarseThre,fineThre);
        end
end

% normalize the given image so as to use uniform thresholding parameters
fff = phi;
fff = fff/max(abs(fff(:)));
fff = fff - mean(fff(:));

%% -------------------------------------------------------------------------
% The synchrosqueezed wave packet transform provids an initial guess
tic;
plotFIG = false;

[radii angles BD TTEng_1st TTEng_2nd masses] = initialGuessTwoStep(fff,subsampleRate,plotFIG,parameters);
% assumes y-coordinate pointing into direction of growing first index
if parameters.BD_type == 1
    [BD agl] = initialGuessOneSector(fff,subsampleRate,plotFIG,parameters2);
    if 0
        figure;imagesc(agl);axis image;colorbar;title('Angles provided by one-sector method');
        figure;imagesc(BD);axis image;title('Boundary provided by one-sector method');
    end
end

angles = mod(angles+pi/18,pi)-pi/18;
% bring wave vectors into order according to angle
[angles,perm] = sort(angles,3);
for k = 1:size(perm,1)
    for l = 1:size(perm,2)
        radii(k,l,:) = radii(k,l,perm(k,l,:));
    end
end

%% Show angles and radii
if 0
    figure;
    subplot(2,3,1);imagesc(angles(:,:,1)*180/pi);axis image;colorbar; title('new angle 1');
    subplot(2,3,2);imagesc(angles(:,:,2)*180/pi);axis image;colorbar; title('new angle 2');
    subplot(2,3,3);imagesc(angles(:,:,3)*180/pi);axis image;colorbar; title('new angle 3');
    subplot(2,3,4);imagesc(radii(:,:,1));axis image;colorbar; title('new radii 1');
    subplot(2,3,5);imagesc(radii(:,:,2));axis image;colorbar; title('new radii 2');
    subplot(2,3,6);imagesc(radii(:,:,3));axis image;colorbar; title('new radii 3');
end

%% -------------------------------------------------------------------------
% find defect region and their connected components
if parameters.BD_type == 2
    temp = smoothImage( sum(masses,3), 4, 1);
    if parameters.isPeriodic == 0
        temp(:,1:2) = 0; temp(:,end-1:end) = 0; temp(1:2,:) = 0; temp(end-1:end,:) = 0;
    end
    binaryIndexCoarse = 1-findDeftArea(temp,energyThre,coarseThre);
    [CCCoarse,colComponentsCoarse] = periodicConnComp( binaryIndexCoarse );
    binaryIndex = 1-findDeftArea(temp,energyThre,fineThre);
else
    pos = find(BD>=BDThre);
    BD(pos) = BDThre;
    if parameters.isPeriodic == 0
        BD(:,1:2) = BDThre; BD(:,end-1:end) = BDThre; BD(1:2,:) = BDThre; BD(end-1:end,:) = BDThre;
    end
    BD = (BDThre-BD)/BDThre*3;
    binaryIndexCoarse = 1-findDeftArea(smoothImage( BD, 4, 1),energyThre,coarseThre);
    [CCCoarse,colComponentsCoarse] = periodicConnComp( binaryIndexCoarse );
    binaryIndex = 1-findDeftArea(smoothImage( BD, 4, 1),energyThre,fineThre);
end

[CC,colComponents] = periodicConnComp( binaryIndex );
[CC2,colComponents2] = periodicConnComp( 1-binaryIndex );

if 0
    pic = figure; imagesc(colComponentsCoarse); axis image; title('coarse defect regions w/ p.b.c.'); colorbar;
    fileName = sprintf('../../results/Part/CoarseDefectEx%d',example);colormap (1-gray);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
    pic = figure; imagesc(colComponents); axis image; title('defect regions w/ p.b.c.'); colorbar;
    fileName = sprintf('../../results/Part/FineDefectEx%d',example);colormap (1-gray);
    set(gca,'xtick',[]);set(gca,'ytick',[]);
end
toc

%% -------------------------------------------------------------------------
tic;
G1 = initGFromSSPeaks( angles, radii, pi/3*[0;1;2], median(radii(:))*ones(3,1));

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

%% -------------------------------------------------------------------------
% show initialization in the whole domain
address = sprintf('results/Part/');
fileName = sprintf('IntWhole');
visualize(phi,G0,curlG0,masses,disVolThre,binaryIndex,0,address,fileName,example, CC);
toc

tic;
%% -------------------------------------------------------------------------
% elastic denoising
G = elasticDenoising( G0, binaryIndex, prescribedCurl, pointGroup, L, numIter );
curlG = findCurl( G, pointGroup, L );

%% -------------------------------------------------------------------------
% show optimized results in the whole domain
address = sprintf('results/Part/');
fileName = sprintf('OptWhole');
visualize(phi,G,curlG,masses,disVolThre,binaryIndex,0,address,fileName,example, CC);
toc
end
