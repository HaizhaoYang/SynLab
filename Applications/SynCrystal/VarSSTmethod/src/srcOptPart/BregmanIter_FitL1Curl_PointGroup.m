function [G,curlG] = BregmanIter_FitL1Curl_PointGroup( waveVecs, masses, stencil, R, G, weights, lambda, nt, GPUflag )
% function [G,curlG] = BregmanIter_FitL1Curl_PointGroup( waveVecs, masses, stencil, R, G, weights, lambda, nt, GPUflag )
% assumes grid spacing of 1
% assumes origin in top left corner
% x-coordinate points right, y-coordinate down
% assumes periodic boundary conditions
% R is [m x n x 2] with rotation angles (expected to correspond to crystal symmetries)
%   for the differences (yp,xp), in that order.
% psi / b are optional arguments. Initialized to curlG and 0, resp., if not given.

  %% set parameters and GPU variables
  format compact;
  maxNumSteps = 5;

  m = size(waveVecs,1);
  if GPUflag
    v_GPU = gpuArray( waveVecs );
    m_GPU = gpuArray( masses );
    G_GPU = gpuArray( G );
    R_GPU = gpuArray( R );
  else
    v_GPU = waveVecs;
    m_GPU = masses;
    G_GPU = G;
    R_GPU = R;
  end

  %% define and initialize psi (the dual variable) and b (the error variable)
  curlG = getCurlR( G_GPU, R_GPU );
  psi = curlG;
  b = 0 * psi;
  
  %% do split Bregman iteration
  for ii = 1:nt
    %% update G (the primal variable)
    fun = @(G) primalEnergyR( G, v_GPU, m_GPU, R_GPU, stencil, psi, b, weights, lambda );
    G_GPU = descendNCG( fun, G_GPU, maxNumSteps, false );
    %% update psi (the dual variable)
    curlG = getCurlR( G_GPU, R_GPU );
    psi = shrink( curlG + b, weights(2) / lambda );
    %% update b (the error variable)
    b = b + curlG - psi;

    %% compute energy value
    [e1,~] = fittingEnergy( G_GPU, v_GPU, m_GPU, stencil );
    [e2,~] = smoothingEnergyR( G_GPU, R_GPU );
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

function curlG = getCurlR( G, R )
  dG2dx  = cat(3,getAdjGComp(G,R,1,2,2),getAdjGComp(G,R,2,2,2)) - G(:,:,:,2);
  dG1dy  = cat(3,getAdjGComp(G,R,1,1,1),getAdjGComp(G,R,2,1,1)) - G(:,:,:,1);
  curlG = -dG1dy+dG2dx;
end

function X = shrink( X, cutoff )
  nrm = repmat(sqrt(X(:,:,1).^2+X(:,:,2).^2),[1 1 2]);
  ix = nrm <= cutoff;
  X(ix) = 0;
  X(~ix) = X(~ix) - X(~ix)./nrm(~ix) * cutoff;
end

function [e,de] = primalEnergyR( G, waveVecs, masses, R, stencil, psi, b, weights, lambda )
  % fitting energy term
  [e,de] = fittingEnergy( G, waveVecs, masses, stencil );
  de = weights(1) * de;
  e = weights(1) * e;
  % smoothing energy term
  [e2,de2] = smoothingEnergyR( G, R );
  de = de + weights(3) * de2;
  e = e + weights(3) * e2;
  % elastic energy term
  [e2,de2] = elasticEnergy( G );
  de = de + weights(4) * de2;
  e = e + weights(4) * e2;
  % Bregman penalty term
  [e2,de2] = bregmanPenaltyR( G, R, psi, b );
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

function [e,de] = smoothingEnergyR( G, R )
  %% compute unscaled gradient of G
  dGdx = ( G - cat(4, cat(3,getAdjGComp(G,R,1,1,2),getAdjGComp(G,R,2,1,2)), ...
                      cat(3,getAdjGComp(G,R,1,2,2),getAdjGComp(G,R,2,2,2))) );
  dGdy = ( G - cat(4, cat(3,getAdjGComp(G,R,1,1,1),getAdjGComp(G,R,2,1,1)), ...
                      cat(3,getAdjGComp(G,R,1,2,1),getAdjGComp(G,R,2,2,1))) );
  %% compute energy
  e = sum(dGdx(:).^2 + dGdy(:).^2);
  %% compute derivative
  de = G;
  de(:,:,1,1) = 2*( dGdx(:,:,1,1)+dGdy(:,:,1,1)-dGdx(:,[end 1:end-1],1,1).*cos(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,1,1).*cos(R([end 1:end-1],:,1))...
                                               -dGdx(:,[end 1:end-1],2,1).*sin(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,2,1).*sin(R([end 1:end-1],:,1)) );
  de(:,:,1,2) = 2*( dGdx(:,:,1,2)+dGdy(:,:,1,2)-dGdx(:,[end 1:end-1],1,2).*cos(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,1,2).*cos(R([end 1:end-1],:,1))...
                                               -dGdx(:,[end 1:end-1],2,2).*sin(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,2,2).*sin(R([end 1:end-1],:,1)) );
  de(:,:,2,1) = 2*( dGdx(:,:,2,1)+dGdy(:,:,2,1)-dGdx(:,[end 1:end-1],2,1).*cos(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,2,1).*cos(R([end 1:end-1],:,1))...
                                               +dGdx(:,[end 1:end-1],1,1).*sin(R(:,[end 1:end-1],2))+dGdy([end 1:end-1],:,1,1).*sin(R([end 1:end-1],:,1)) );
  de(:,:,2,2) = 2*( dGdx(:,:,2,2)+dGdy(:,:,2,2)-dGdx(:,[end 1:end-1],2,2).*cos(R(:,[end 1:end-1],2))-dGdy([end 1:end-1],:,2,2).*cos(R([end 1:end-1],:,1))...
                                               +dGdx(:,[end 1:end-1],1,2).*sin(R(:,[end 1:end-1],2))+dGdy([end 1:end-1],:,1,2).*sin(R([end 1:end-1],:,1)) );
end

function [e,de] = bregmanPenaltyR( G, R, psi, b )
  %% compute Bregman term
  curlG = getCurlR( G, R );
  bregmanSum = curlG - psi + b;
  %% compute energy
  e = 1 / 2 * sum( bregmanSum(:).^2 );
  %% compute derivative
  de = G;
  de(:,:,1,1) =  bregmanSum(:,:,1) - bregmanSum([end 1:end-1],:,1) .* cos(R([end 1:end-1],:,1)) - bregmanSum([end 1:end-1],:,2) .* sin(R([end 1:end-1],:,1));
  de(:,:,1,2) = -bregmanSum(:,:,1) + bregmanSum(:,[end 1:end-1],1) .* cos(R(:,[end 1:end-1],2)) + bregmanSum(:,[end 1:end-1],2) .* sin(R(:,[end 1:end-1],2));
  de(:,:,2,1) =  bregmanSum(:,:,2) + bregmanSum([end 1:end-1],:,1) .* sin(R([end 1:end-1],:,1)) - bregmanSum([end 1:end-1],:,2) .* cos(R([end 1:end-1],:,1));
  de(:,:,2,2) = -bregmanSum(:,:,2) - bregmanSum(:,[end 1:end-1],1) .* sin(R(:,[end 1:end-1],2)) + bregmanSum(:,[end 1:end-1],2) .* cos(R(:,[end 1:end-1],2));
end

function RG = getAdjGComp( G, R, ii, jj, kk )
% function RG = getAdjGComp( G, R, ii, jj, kk )
% Computes G(:,:,ii,jj) with the shift specified by kk,
% (yp,xp), resp.

  if(kk == 1)
    if(ii == 1)
      if(jj == 1)
        RG = cos(R(:,:,kk)).*G([2:end 1],:,1,1)-sin(R(:,:,kk)).*G([2:end 1],:,2,1);
      else % jj == 2
        RG = cos(R(:,:,kk)).*G([2:end 1],:,1,2)-sin(R(:,:,kk)).*G([2:end 1],:,2,2);
      end
    else % ii == 2
      if(jj == 1)
        RG = sin(R(:,:,kk)).*G([2:end 1],:,1,1)+cos(R(:,:,kk)).*G([2:end 1],:,2,1);
      else % jj == 2;
        RG = sin(R(:,:,kk)).*G([2:end 1],:,1,2)+cos(R(:,:,kk)).*G([2:end 1],:,2,2);
      end
    end
  else % (kk == 2)
    if(ii == 1)
      if(jj == 1)
        RG = cos(R(:,:,kk)).*G(:,[2:end 1],1,1)-sin(R(:,:,kk)).*G(:,[2:end 1],2,1);
      else % jj == 2
        RG = cos(R(:,:,kk)).*G(:,[2:end 1],1,2)-sin(R(:,:,kk)).*G(:,[2:end 1],2,2);
      end
    else % ii == 2
      if(jj == 1)
        RG = sin(R(:,:,kk)).*G(:,[2:end 1],1,1)+cos(R(:,:,kk)).*G(:,[2:end 1],2,1);
      else % jj == 2;
        RG = sin(R(:,:,kk)).*G(:,[2:end 1],1,2)+cos(R(:,:,kk)).*G(:,[2:end 1],2,2);
      end
    end
  end
  
end
