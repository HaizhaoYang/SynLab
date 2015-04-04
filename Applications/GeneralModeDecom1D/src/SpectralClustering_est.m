%--------------------------------------------------------------------------
% This function takes a NxN matrix CMat as adjacency of a graph and 
% computes the segmentation of data from spectral clustering. It estimates
% the number of subspaces using eigengap heuristic
% CMat: NxN adjacency matrix
% K: number of largest coefficients to choose from each column of CMat
% Grps: [grp1,grp2,grp3] for three different forms of Spectral Clustering
% SingVals: [SV1,SV2,SV3] singular values for three different forms of SC
% LapKernel(:,:,i): n last columns of kernel of laplacian to apply KMeans
%--------------------------------------------------------------------------
% By Mahdi Soltanolkotabi and Emmanuel Candes
% Modified version of SpectralClustering written by Ehsan Elhamifar
%--------------------------------------------------------------------------

function [group , svals, Lap , n ] = SpectralClustering_est(CKSym, option )


N = size(CKSym,1);
MAXiter = 1000; % Maximum iteration for KMeans Algorithm
REPlic = 100;   % Replication for KMeans Algorithm

switch option
    
    case 1
        % Method 1: Unnormalized Method
        DKU = diag( sum(CKSym) );
        Lap  = DKU - CKSym;
        [uKU,sKU,vKU] = svd(Lap );
        
        svals = diag(sKU);
        [ min_val , ind_min ] = min( diff( svals(1:end-1) ) ) ;  
        n = size(CKSym , 1 ) - ind_min ;
        
        f = size(vKU,2);
        kerKU = vKU(:,f-n+1:f);
        
        group = kmeans(kerKU,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
        
    case 2
        % Method 2: Random Walk Method
        DKN=( diag( sum(CKSym) ) )^(-1);
        Lap = speye(N) - DKN * CKSym;
        [uKN,sKN,vKN] = svd( Lap );
        svals = diag(sKN);
        [ min_val , ind_min ] = min( diff( svals(1:end-1) ) ) ;  
        n = size(CKSym , 1 ) - ind_min ;
        
        f = size(vKN,2);
        kerKN = vKN(:,f-n+1:f);
        

        group = kmeans(kerKN,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
        
    case 3
        % Method 3: Normalized Symmetric
        DKS = ( diag( sum(CKSym) ) )^(-1/2);
        Lap  = speye(N) - DKS * CKSym * DKS;
        [uKS,sKS,vKS] = svd(Lap );
        svals = diag(sKS);
       
        [ min_val , ind_min ] = min( diff( svals(1:end-1) ) ) ;  
        n = size(CKSym , 1 ) - ind_min ;
        
        f = size(vKS,2);
        kerKS = vKS(:,f-n+1:f);
        for i = 1:N
            kerKS(i,:) = kerKS(i,:) ./ norm(kerKS(i,:));
        end
        
        group = kmeans(kerKS,n,'start','sample','maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
end 