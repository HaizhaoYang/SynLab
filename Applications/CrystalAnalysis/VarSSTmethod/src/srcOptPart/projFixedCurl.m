function [projG,cache] = projFixedCurl( G, R, L, Curl, cache )
% G     - 2D array of matrices, M x N x 2 x 2
% R     - list of all point group elements (i.e. rotation matrices), 2 x 2 x H
% L     - list of pixel pairs across which we use a nontrivial point group element, K x 4,
%         row entries are row and column index of first pixel, right or up difference, index of point group element
% Curl  - given fixed curl for projection (empty counts as zero)
% cache - contains auxiliary data that has to be precomputed only once for fixed R,L; set cache.valid=false for fresh computation
% projects G orthogonally onto the constrained set findCurl(G)=Curl
curlG = findCurl(G,R,L);
if ~isempty( Curl )
    curlG = curlG - Curl;
end
[V,cache] = LaplInv(curlG,R,L,cache);
projG = G - curlAdj(V,R,L);
end

function test
% create test data
M = 128; N = M;
G = randn(M,N,2,2);
rot = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
R = eye(2);
for j = 1:5
    R = cat(3,R,R(:,:,j)*rot);
end
K = 20;
L = [ceil(M*rand(K,1)),ceil(N*rand(K,1)),round(rand(K,1)),ceil(6*rand(K,1))];
cache.valid = false;

% check consistency of curl and curlAdj
curlG = findCurl(G,R,L);
V = rand(M,N,2);
curlAdjV = curlAdj(V,R,L);
sum(G(:)'*curlAdjV(:))-sum(V(:)'*curlG(:))

% check consistency of curl, curlAdj, and LaplInv
[V,cache] = LaplInv(findCurl(curlAdj(V,R,L),R,L),R,L,cache);
NegV = -LaplInv(findCurl(curlAdj(V,R,L),R,L),R,L,cache);
norm(NegV(:)+V(:),2)

% check projection onto zero curl
projG = G - curlAdj(LaplInv(curlG,R,L,cache),R,L);
curlProjG = findCurl(projG,R,L);
norm(curlProjG(:),2)
end

function curlAdjV = curlAdj( V, R, L )
% V - 2D array of vectors, M x N x 2
% R - list of all point group elements (i.e. rotation matrices), 2 x 2 x H
% L - list of pixel pairs across which we use a nontrivial point group element, K x 4,
%     row entries are row and column index of first pixel, right or up difference, index of point group element
% computes adjoint of curl

[M,N,~] = size(V);

% compute normal curl adjoint
curlAdjV = cat(4,V-V([2:end 1],:,:),V(:,[end 1:end-1],:)-V);

% perturb by point group operator
RminI = R;
RminI(1,1,:) = RminI(1,1,:)-1;
RminI(2,2,:) = RminI(2,2,:)-1;
for k = 1:size(L,1)
    aux = num2cell(L(k,:));
    [m,n,t,s] = deal(aux{:});
    if t
        curlAdjV(m,mod(n,N)+1,:,2) = curlAdjV(m,mod(n,N)+1,:,2) + shiftdim( RminI(:,:,s)' * squeeze(V(m,n,:)), -2 );
    else
        curlAdjV(mod(m-2,M)+1,n,:,1) = curlAdjV(mod(m-2,M)+1,n,:,1) - shiftdim( RminI(:,:,s)' * squeeze(V(m,n,:)), -2 );
    end
end

end

function [LaplInvV,cache] = LaplInv( V, R, L, cache )
% V - 2D array of vectors, M x N x 2
% R - list of all point group elements (i.e. rotation matrices), 2 x 2 x H
% L - list of pixel pairs across which we use a nontrivial point group element, K x 4,
%     row entries are row and column index of first pixel, right or up difference, index of point group element
% cache - contains auxiliary data that has to be precomputed only once for fixed R,L
% solves -Laplace LaplInvV = V (i.e. returns one possible solution if solvable)

[M,N,~] = size(V);
if ~cache.valid
    cache.valid = true;
    
    [kx,ky] = meshgrid(0:N-1,0:M-1);
    ddx = cos(2*pi/N*kx)+i*sin(2*pi/N*kx) - 1; % fft of [-1 0 ... 0 1;0 ...] (note: convolution yields forward difference)
    ddy = cos(2*pi/M*ky)-i*sin(2*pi/M*ky) - 1; % fft of [-1 0 ...;1 0 ...;0 ...]
    cache.Laplace = abs(ddx).^2 + abs(ddy).^2; cache.Laplace(1,1) = Inf; % dividing by Laplace now always kills the constant term; actually is -Laplace
    
    % if not yet available, produce the point group operator J_\Delta in sparse matrix form and identify the support of its range
    RminI = R;
    RminI(1,1,:) = RminI(1,1,:)-1;
    RminI(2,2,:) = RminI(2,2,:)-1;
    rows = zeros(2,2,size(L,1)); cols = rows; vals = rows;
    for k = 1:size(L,1)
        aux = num2cell(L(k,:));
        [m,n,t,s] = deal(aux{:});
        indices = sub2ind([M,N],[m mod(m-2,M)+1 m],[n n mod(n,N)+1]);
        rows(:,:,k) = [indices(1) indices(1);M*N+indices(1) M*N+indices(1)];
        cols(:,:,k) = [indices(2+t) M*N+indices(2+t);indices(2+t) M*N+indices(2+t)];
        vals(:,:,k) = RminI(:,:,s);
    end
    cache.J_Delta = sparse([rows(:);cols(:)],[cols(:);rows(:)],[vals(:);vals(:)],2*M*N,2*M*N); % this is the operator applied to V(:)
    
    % if not yet available, identify a basis s_1,s_2,... of range J_\Delta; each s_i only contains a single 1;
    % K x 4, where ith row is linear, row, column, and 3rd dim index of nonzero element of s_i
    S = find(sum(cache.J_Delta))';
    [I,J,H] = ind2sub([M N 2],S);
    cache.S = [S I J H];
    cache.K = length(I);
    
    % if not yet available, set up matrix of reduced system and compute QRP decomposition
    e1 = zeros(M,N); e1(1,1) = 1;
    LaplInvE1 = real(ifft2(fft2(e1)./cache.Laplace));
    A = zeros(cache.K+2,cache.K+2);
    for ii = 1:cache.K % ith column of A
        [rowJ,~,valJ] = find(cache.J_Delta(:,cache.S(ii,1)));
        [I,J,H] = ind2sub([M N 2],rowJ);
        for jj = 1:cache.K % jth row of A
            if 0
                for k = 1:size(rowJ,1)
                    if H(k) == cache.S(jj,4)
                        A(jj,ii) = A(jj,ii) + valJ(k) * LaplInvE1(mod(I(k)-cache.S(jj,2),M)+1,mod(J(k)-cache.S(jj,3),N)+1);
                    end
                end
            else
                pos = find(H==cache.S(jj,4));
                A(jj,ii) = sum(valJ(pos).*LaplInvE1(sub2ind([M N],mod(I(pos)-cache.S(jj,2),M)+1,mod(J(pos)-cache.S(jj,3),N)+1)) );
            end
        end
        A(ii,ii) = A(ii,ii) - 1;
        A(end-1,ii) = sum((2-H).*valJ)/M/N;
        A(end,  ii) = sum((H-1).*valJ)/M/N;
        % kth row of A
        A(ii,end+cache.S(ii,4)-2) = 1/M/N;
    end
    [cache.Q,cache.U,Pt] = qr(A,'vector');
    cache.P(Pt) = 1:length(Pt);
    cache.ind = find(~diag(cache.U),1,'first');
    if ~isempty(cache.ind)
        cache.U(cache.ind:end,cache.ind:end) = eye(cache.K+3-cache.ind);
    end
end

% set up rhs of reduced system
VFFT = fft2(V);
VFFT(:,:,1) = VFFT(:,:,1) ./ cache.Laplace;
VFFT(:,:,2) = VFFT(:,:,2) ./ cache.Laplace;
LaplInvV = real(ifft2(VFFT));
rhs = zeros(cache.K+2,1);
for k = 1:cache.K
    rhs(k) = -LaplInvV(cache.S(k,2),cache.S(k,3),cache.S(k,4));
end
rhs(end-1:end) = -mean(mean(V(:,:,:),1),2);

% assemble solution
rhs = cache.Q'*rhs;
if ~isempty(cache.ind)
    rhs(cache.ind:end) = 0;
end
lambda = cache.U\rhs;
lambda = lambda(cache.P);
LaplInvV(:,:,1) = LaplInvV(:,:,1) + lambda(end-1)/M/N;
LaplInvV(:,:,2) = LaplInvV(:,:,2) + lambda(end)/M/N;
weights = cache.J_Delta(cache.S(:,1),cache.S(:,1)) * lambda(1:end-2);
aux = zeros(M,N,2);
aux(cache.S(:,1)) = weights;
LaplInvV = LaplInvV + real(ifft2(fft2(aux)./repmat(cache.Laplace,[1 1 2])));
end