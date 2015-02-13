function [noise_data] = RFI_MakeDataAlphaStable(alpha, disper, delta, N, beta, NDim)
% [noise_data] = RFI_MakeDataAlphaStable(alpha, disper, delta, N, beta, NDim);
% This function generates a vector of synthetic data samples that follow the symmetric alpha
% stable model. When NDim > 1, it generates isotropic symmetric alpha stable model in NDim dimensions.
% Inputs:   alpha  - Characteristic Exponent
%           disper - Dispersion Parameter
%           delta  - Localization Parameter
%           N      - No of samples to be generated
% Optional Inputs:
%           beta   - Skewness Parameter (default = 0)
%           NDim   - Dimension of data begin generated (default = 1)
% Outputs:  noise_data - NDim x N size matrix of noise data generated based on
%                        symmetric alpha stable model. When NDim > 1, it is
%                        restricted to be an isotropic symmetric alpha stable
%                        model in NDim dimensions.
%
% Reference:
% [1] Ercan Kuruoglu, "Signal Processing in alpha-stable noise environments
%     - A least lp-NORM approach", PhD Thesis, University of Cambridge, Nov
%     1998.
% [2] Panagiotis Tsakalides, "Array Signal Processing with Alpha Stable
%     Distributions", PhD Thesis, University of Southern California, 1995
%
% Copyright (c) The University of Texas
% Please see the file Copyright.txt that came with this release for details
% Programmers: Kapil Gulati   (gulati@ece.utexas.edu)
%
% Revisions:
% Date      Version     Owner     Comments
% -------------------------------------------------------------------------
% 09/22/07  1.0         gulati    Initial Creation
% 11/14/07  1.1 Beta    gulati    Fixed the fact that S((alpha,disper,0)
%                                 = disper ^ (1/alpha)* S(alpha, 1, 0) in
%                                 the implementation. Now sync with [1].
% 12/05/09  1.4 Beta    gulati    Added skewness parameter beta
% 12/11/10  1.6 Beta    gulati    Extended function to support generation
%                                 of multi-dimensional isotropic symmetric
%                                 alpha stable random vector.

%Assigned default value
if nargin < 5
    beta = 0;
    NDim = 1;
end;

if nargin < 6
    NDim = 1;
end

% Check for input parameters
if (NDim > 1 && (beta ~= 0 || delta ~= 0))
    error('For NDim > 1, this function supports generation of only isotropic symmetric alpha stable random vector (beta = 0)');
end

if (NDim > 1)
    % Generate samples for isotropic symmetric alpha stable random vector
    % Using method:
    % X = A^(1/2)(G_1, G_2, ..., G_NDim) where X is the NDim dimensional
    % isotropic alpha stable random vector, A is an alpha/2-positive stable
    % random  variable with dispersion 2^(alpha/2)*disper*cos(pi*alpha/4),
    % and G_1, ..., G_NDim are independent Gaussian random variables with
    % zero mean and unit variance.
    A           =  RFI_MakeDataAlphaStable(alpha/2, 2^(alpha/2)*disper*cos(pi*alpha/4), 0, N, 1, 1);
    noise_data  =  repmat(A.^(1/2), NDim, 1).* (mvnrnd(zeros(N, NDim), eye(NDim, NDim)).');
    
    
else
    % NDim = 1
    
    % Refer to [1],[2] for details on this method of data generation.
    if (alpha == 2)         % Gaussian
        noise_data = sqrt(2)*randn(1,N);
    elseif (alpha == 1)     % Cauchy
        u = pi*(rand(1,N)-.5*ones(1,N));
        noise_data = 2/pi * ((pi/2 + beta.*u).*tan(u) - beta* log ((pi/2 * (-log (rand(1,N))) .* cos(u))./(pi/2 + beta*u)));
    else
        u = pi*(rand(1,N)-.5*ones(1,N));
        u0 = -pi/2*beta*(1-abs(1-alpha))/alpha;
        noise_data = (cos(atan(beta * (tan(pi*alpha/2)))))^(-1/alpha)*(sin (alpha.*(u - u0))) ./ ((cos(u)).^(1/alpha)) .* ((cos((1-alpha).*u + alpha.*u0)) ./ (-log (rand(1,N)))).^((1-alpha)/alpha);
    end
    %
    
    if (alpha == 1)
        noise_data = (disper)^(1/alpha) .* noise_data + 2/pi * beta  * disper^(1/alpha) * 1/alpha * log (disper) + delta;      % Including dispersion and localization
    else
        noise_data = (disper)^(1/alpha) .* noise_data + delta;      % Including dispersion and localization
    end
end

return