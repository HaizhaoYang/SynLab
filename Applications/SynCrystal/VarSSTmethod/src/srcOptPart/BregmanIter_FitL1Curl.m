function [G,curlG] = BregmanIter_FitL1Curl( waveVecs, masses, stencil, G, weights, lambda, nt, GPUflag )
% function [G,curlG] = l1CurlBregman( waveVecs, masses, stencil, G, weights, lambda, nt1, GPUflag )
% assumes grid spacing of 1
% assumes origin in top left corner
% x-coordinate points right, y-coordinate down
% assumes periodic boundary conditions

  %% set parameters and GPU variables
  format compact;
  maxNumSteps = 5;

  m = size(waveVecs,1);
  if GPUflag
    v_GPU = gpuArray( waveVecs );
    m_GPU = gpuArray( masses );
    G_GPU = gpuArray( G );
  else
    v_GPU = waveVecs;
    m_GPU = masses;
    G_GPU = G;
  end

  %% define and initialize psi (the dual variable) and b (the error variable)
  curlG = getCurl( G_GPU );
  psi = curlG;
  b = 0 * psi;

  %% do split Bregman iteration
  for ii = 1:nt
    %% update G (the primal variable)
    fun = @(G) primalEnergy( G, v_GPU, m_GPU, stencil, psi, b, weights, lambda );
    G_GPU = descendNCG( fun, G_GPU, maxNumSteps, false );
    %% update psi (the dual variable)
    curlG = getCurl( G_GPU );
    psi = shrink( curlG + b, weights(2) / lambda );
    %% update b (the error variable)
    b = b + curlG - psi;

    %% compute energy value
    [e1,~] = fittingEnergy( G_GPU, v_GPU, m_GPU, stencil );
    [e2,~] = smoothingEnergy( G_GPU );
    [e3,~] = elasticEnergy( G_GPU );
    fittE = weights(1) * e1;
    l1cuE = weights(2) * sum( sum( sqrt( sum( curlG.^2, 3 ) ) ) );
    smthE = weights(3) * e2;
    elasE = weights(4) * e3;
    if(GPUflag)
      fittE = gather( fittE );
      l1cuE = gather( l1cuE );
      smthE = gather( smthE );
      elasE = gather( elasE );
    end
    fprintf( 'step %d, total energy %f, fitting energy %f, L1-norm of curl %f, smoothing energy %f, elastic energy %f\n', ii, fittE + l1cuE + smthE + elasE, fittE, l1cuE, smthE, elasE )
  end
  
  %% return result
  if GPUflag
    G     = gather( G_GPU );
    curlG = gather( curlG );
  else
    G = G_GPU;
  end
end

function curlG = getCurl( G )
  dG2dx  = G(:,[2:end 1],:,2) - G(:,:,:,2);
  dG1dy  = G([2:end 1],:,:,1) - G(:,:,:,1);
  %dG2dx  = circshift( G(:,:,:,2), [0 -1 0] ) - G(:,:,:,2);
  %dG1dy  = circshift( G(:,:,:,1), [-1 0 0] ) - G(:,:,:,1);
  curlG = dG2dx - dG1dy;
end

function X = shrink( X, cutoff )
  nrm = repmat(sqrt(X(:,:,1).^2+X(:,:,2).^2),[1 1 2]);
  ix = nrm <= cutoff;
  X(ix) = 0;
  X(~ix) = X(~ix) - X(~ix)./nrm(~ix) * cutoff;
end

function [e,de] = primalEnergy( G, waveVecs, masses, stencil, psi, b, weights, lambda )
  % fitting energy term
  [e,de] = fittingEnergy( G, waveVecs, masses, stencil );
  de = weights(1) * de;
  e = weights(1) * e;
  % smoothing energy term
  [e2,de2] = smoothingEnergy( G );
  de = de + weights(3) * de2;
  e = e + weights(3) * e2;
  % elastic energy term
  [e2,de2] = elasticEnergy( G );
  de = de + weights(4) * de2;
  e = e + weights(4) * e2;
  % Bregman penalty term
  [e2,de2] = bregmanPenalty( G, psi, b );
  de = de + lambda * de2;
  e = e + lambda * e2;
end

function [e,de] = fittingEnergy( G, waveVecs, masses, stencil )
  e = 0;
  de = 0 * G;
  %% run over all stencil points
  for k = 1:size(stencil,1)
    %% compute G^T stencil
    v1 = G(:,:,1,1) * stencil(k,1) + G(:,:,2,1) * stencil(k,2);
    v2 = G(:,:,1,2) * stencil(k,1) + G(:,:,2,2) * stencil(k,2);
    %% compute mismatch, energy density, and gradient
    mismatch = cat(4,v1,v2) - waveVecs(:,:,k,:);
    e = e + sum(sum(masses(:,:,k).*sum(mismatch.^2,4)));
    de(:,:,1,1) = de(:,:,1,1) + 2 * stencil(k,1) * masses(:,:,k).*mismatch(:,:,:,1);
    de(:,:,2,1) = de(:,:,2,1) + 2 * stencil(k,2) * masses(:,:,k).*mismatch(:,:,:,1);
    de(:,:,1,2) = de(:,:,1,2) + 2 * stencil(k,1) * masses(:,:,k).*mismatch(:,:,:,2);
    de(:,:,2,2) = de(:,:,2,2) + 2 * stencil(k,2) * masses(:,:,k).*mismatch(:,:,:,2);
  end
end

function [e,de] = elasticEnergy( G )
% implements a nonlinear elastic energy with cutoff;
% here: arctan(W(G)) with W(G) = mu/2|G|^2+lambda/4*(det G)^2-(mu+lambda/2)log(det G)-mu-lambda/4
  mu = 1;
  lambda = 1;
  detG = G(:,:,1,1) .* G(:,:,2,2) - G(:,:,1,2) .* G(:,:,2,1);
  if min(detG(:)) <= 0
    e = Inf;
    de = Inf;
    return;
  end
  enDens = - mu - lambda/4 + mu/2*sum(sum(G.^2,4),3) + lambda/4*detG.^2 - (mu+lambda/2)*log(detG);
  e = sum(sum( atan(enDens) ));
  enDens = 1 ./ ( 1 + enDens.^2 ); % reuse variable enDens
  detG = lambda/2*detG - (mu+lambda/2)./detG; % reuse variable detG
  de = G;
  de(:,:,1,1) = enDens .* ( mu * G(:,:,1,1) + detG .* G(:,:,2,2) );
  de(:,:,2,1) = enDens .* ( mu * G(:,:,2,1) - detG .* G(:,:,1,2) );
  de(:,:,1,2) = enDens .* ( mu * G(:,:,1,2) - detG .* G(:,:,2,1) );
  de(:,:,2,2) = enDens .* ( mu * G(:,:,2,2) + detG .* G(:,:,1,1) );
end

function [e,de] = smoothingEnergy( G )
  %% compute unscaled gradient of G
  dGdx  = ( G(:,[2:end 1],:,:) - G );
  dGdy  = ( G([2:end 1],:,:,:) - G );
  %dGdx  = ( circshift( G, [0 -1 0 0] ) - G );
  %dGdy  = ( circshift( G, [-1 0 0 0] ) - G );
  %% compute energy
  e = sum( sum( sum( sum( dGdx.^2 + dGdy.^2 ) ) ) );
  %% compute derivative
  de = dGdx(:,[end,1:end-1],:,:) + dGdy([end,1:end-1],:,:,:) - dGdx - dGdy;
  %de = circshift( dGdx, [0 1 0 0] ) + circshift( dGdy, [1 0 0 0] ) - dGdx - dGdy;
  de = 2 * de;
end

function [e,de] = bregmanPenalty( G, psi, b )
  %% compute Bregman term
  curlG = getCurl( G );
  bregmanSum = curlG - psi + b;
  %% compute energy
  e = 1 / 2 * sum( sum( sum( bregmanSum.^2 ) ) );
  %% compute derivative
  de = G;
  de(:,:,:,2) = bregmanSum(:,[end,1:end-1],:) - bregmanSum;
  de(:,:,:,1) = bregmanSum - bregmanSum([end,1:end-1],:,:);
  %de(:,:,:,2) = circshift( bregmanSum, [0 1 0] ) - bregmanSum;
  %de(:,:,:,1) = bregmanSum - circshift( bregmanSum, [1 0 0] );
end
