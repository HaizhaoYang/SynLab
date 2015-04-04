function [x, cost] = TVD_mm(y, lam, Nit)
% [x, cost] = TVD_mm(y, lam, Nit)
% Total variation denoising using majorization-minimization
% and a sparse linear system solver.
%
% INPUT
%   y - noisy signal
%   lam - regularization parameter
%   Nit - number of iterations
%
% OUTPUT
%   x - denoised signal
%   cost - cost function history
%
% Reference
% 'On total-variation denoising: A new majorization-minimization
% algorithm and an experimental comparison with wavalet denoising.'
% M. Figueiredo, J. Bioucas-Dias, J. P. Oliveira, and R. D. Nowak.
% Proc. IEEE Int. Conf. Image Processing, 2006.

% Ivan Selesnick
% selesi@poly.edu
% 2011

y = y(:);                                               % ensure column vector
cost = zeros(1, Nit);                                   % cost function history
N = length(y);

e = ones(N-1, 1);
DDT = spdiags([-e 2*e -e], [-1 0 1], N-1, N-1);         % D*D' (sparse matrix)
D = @(x) diff(x);                                       % D (operator)
DT = @(x) [-x(1); -diff(x); x(end)];                    % D'

x = y;                                                  % initialization
Dx = D(x);

for k = 1:Nit
    F = 2/lam * spdiags(abs(Dx), 0, N-1, N-1) + DDT;    % F : sparse matrix
    % F = 2/lam * diag(abs(D(x))) + DDT;                % not a sparse matrix

    x = y - DT(F\D(y));                                 % Solve sparse linear system
    Dx = D(x);
    
    cost(k) = sum(abs(x-y).^2) + lam * sum(abs(Dx));    % keep cost function history
end
