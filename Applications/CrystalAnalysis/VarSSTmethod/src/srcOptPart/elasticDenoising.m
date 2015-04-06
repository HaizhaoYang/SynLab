function G = elasticDenoising( G0, mask, B, R, L, numIter )
% Minimizes the functional
% E[G]=\int_{\Omega\setminus D} (G-G_0)^2 + hyperelastic energy(G) dx
% under the constraint findCurl( G = \sum_{k=1}^K b_k \delta_{x_k}.
% Assumes periodic image pixel spacing 1.
% G0   - G_0 (m x n x 2 x 2)
% mask - characteristic function of \Omega\setminus D (m x n)
% B    - prescribed findCurl( image
% R    - elements of crystal point group (2 x 2 x M; each R(:,:,i) is a rotation matrix)
% L    - list of pixel pairs across which, when we compute a finite difference, we multiply G at the second pixel with a nontrivial point group element,
%        N x 4, row entries are row and column index of first pixel, right or up difference, index of point group element
  
  [m,n,~,~] = size(G0);

  % produce feasible initial guess
  cache.valid = false;
  [G,cache] = projFixedCurl( G0, R, L, B, cache );

  % perform optimization
  proj = @(G) projFixedCurl( G, R, L, [], cache );
  fun = @(G) denoisingEnergy( G, G0, mask, false );
  G = descendNCGConstr( fun, G, proj, numIter, 1 );
end

function [e,de] = denoisingEnergy( G, G0, mask, vecReshaped )
% implements squared L2-distance from G0 plus elastic energy
  if vecReshaped
    G = reshape(G,size(G0));
  end
  [e,de] = elasticEnergy( G, mask );
  e = e + sum(sum( mask .* sum(sum( (G-G0).^2, 4 ), 3 ) ));
  for i = 1:2
    for j = 1:2
      de(:,:,i,j) = de(:,:,i,j) + 2*mask .* (G(:,:,i,j)-G0(:,:,i,j));
    end
  end
  if vecReshaped
    de = de(:);
  end
end

function [e,de] = elasticEnergy( G, mask )
% implements a nonlinear elastic energy
% here: W(G) = mu/2|G|^2+lambda/4*(det G)^2-(mu+lambda/2)log(det G)-mu-lambda/4
% or, to have something defined for neg. determinant: W(G) = mu/2|G|^2+lambda/4*(det G)^2+(mu+lambda/2)(det G-1)(det G-3)/2-mu-lambda/4
  mu = 1;
  lambda = 1;
  detG = G(:,:,1,1) .* G(:,:,2,2) - G(:,:,1,2) .* G(:,:,2,1);
  if 0 & ( min(detG(:)) <= 0 )
    e = Inf;
    de = Inf;
    return;
  end
  %enDens = - mu - lambda/4 + mu/2*sum(sum(G.^2,4),3) + lambda/4*detG.^2 - (mu+lambda/2)*log(detG);
  enDens = - mu - lambda/4 + mu/2*sum(sum(G.^2,4),3) + lambda/4*detG.^2 + (mu+lambda/2)/2*(detG-1).*(detG-3);
  e = sum(sum(mask.*enDens));
  %detG = lambda/2*detG - (mu+lambda/2)./detG; % reuse variable detG
  detG = lambda/2*detG + (mu+lambda/2)*(detG-2); % reuse variable detG
  de = G;
  de(:,:,1,1) = mask .* ( mu * G(:,:,1,1) + detG .* G(:,:,2,2) );
  de(:,:,2,1) = mask .* ( mu * G(:,:,2,1) - detG .* G(:,:,1,2) );
  de(:,:,1,2) = mask .* ( mu * G(:,:,1,2) - detG .* G(:,:,2,1) );
  de(:,:,2,2) = mask .* ( mu * G(:,:,2,2) + detG .* G(:,:,1,1) );
end

