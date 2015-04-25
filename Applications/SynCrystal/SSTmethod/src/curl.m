function curlG = curl( G, R, L )
% G - 2D array of matrices, M x N x 2 x 2
% R - list of all point group elements (i.e. rotation matrices), 2 x 2 x H
% L - list of pixel pairs across which we use a nontrivial point group element, K x 4,
%     row entries are row and column index of first pixel, right or up difference, index of point group element
% computes curl of the rows of a matrix field, where x- and y-direction are in direction of increasing second and decreasing first index, respectively

  % compute normal curl
  curlG = G(:,[2:end 1],:,2) - G(:,:,:,2) + G(:,:,:,1) - G([end 1:end-1],:,:,1);

  % perturb by point group operator
  if nargin >= 3
    RminI = R;
    RminI(1,1,:) = RminI(1,1,:)-1;
    RminI(2,2,:) = RminI(2,2,:)-1;
    for k = 1:size(L,1)
      aux = num2cell(L(k,:));
      [m,n,t,s] = deal(aux{:});
      if t
        curlG(m,n,:) = curlG(m,n,:) + shiftdim( RminI(:,:,s) * squeeze(G(m,mod(n,end)+1,:,2)), -2 );
      else
        curlG(m,n,:) = curlG(m,n,:) - shiftdim( RminI(:,:,s) * squeeze(G(mod(m-2,end)+1,n,:,1)), -2 );
      end
    end
  end
  
end
