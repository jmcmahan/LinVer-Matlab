function [t, w] = gauss_quadrature(measure, order, mparam)
% [t, w] = gauss_quadrature(measure, order) - Compute the Gaussian
% quadrature nodes and weights for a selected integration order and
% measure for integration. For the uniform and beta measures, 
% integration is over [-1, 1], for the gaussian measure, it is 
% -infty, infty. 
%
% measure = 'uniform', 'beta', 'gaussian', 'gamma', a string selecting 
% which  probability measure the integral is weighted by. 
%
% order = order of the integration (number of nodes)
%
% mparam = supplemental parameters for the measures. Use for the
% following measures:
% beta (note this is the (1-x)^alpha * (1+x)^beta notation)
% **********************************************************
% mparam.alpha = alpha parameter
% mparam.beta = beta parameter
% 
% gamma (this is the x^(alpha-1)*exp(-x/beta) version)
% **********************************************************
% mparam.alpha = alpha parameter
% mparam.beta = beta parameter
%
% Returns:
% t - nodes to evalate the integrand on
% w - weights for the quadrature rule

N = order; 

if strcmp(measure, 'gaussian')
    % Recurrence relation for the Hermite polynomials
    a = ones(N, 1);
    b = zeros(N, 1);
    c = (0:(N-1))';
elseif strcmp(measure, 'uniform')
    % Recurrence relation for the Legendre polynomials
    n = N-1; 
    a = (2*(0:n) +1)'./((0:n)+1)';
    b = zeros(N,1);
    c = (0:n)'./((0:n)+1)';
elseif strcmp(measure, 'beta')
    % Recurrence relation for the Jacobi polynomials
    if nargin < 3
        al = 1; bt = 1; 
    else
        al = mparam.alpha;
        bt = mparam.beta;
    end

    n = (0:(N-1)); 
    a = 0.5 * (2*n+al+bt+1) .* (2*n+al+bt+2) ...
        ./ ((n+1) .* (n+al+bt+1));
    b = 0.5 * (al^2-bt^2) .* (2*n+al+bt+1) ...
        ./ ((n+1) .* (n+al+bt+1) .* (2*n+al+bt));
    c =       (n+al).*(n+bt).*(2*n+al+bt+2) ...
        ./ ((n+1) .* (n+al+bt+1) .* (2*n+al+bt));

    % Correction to first term
    a(1) = 0.5 * (al + bt + 2);
    b(1) = 0.5 * (al - bt); 
elseif strcmp(measure, 'gamma')
    % Recurrence relation for the Laguerre polynomials
    if nargin < 3
        al = 2; bt = 1;
    else
        al = mparam.alpha;
        bt = mparam.beta;
    end
    n = (0:(N-1));
    a = -1./(n+1);
    b = (2*n + al)./(n+1);
    c = n./(al -1 + n);
end


T = -diag(b./a) + diag(1./a(1:N-1), 1) + diag(c(2:N)./a(2:N),-1);


if sum(sum(T ~= T'))
    alpha = -b ./ a;
    beta = sqrt(c(2:end)./(a(1:(N-1)).*a(2:end)));
    T = diag(alpha) + diag(beta, 1) + diag(beta, -1); 
end

[V, t] = eig(T); 

% Nodes / abscissa
t = diag(t);
% Weights
w = ((V(1,:)).^2)';


if strcmp(measure, 'gamma')
    % Scale the nodes in this special case
    t = bt*t; 
end

end


function scratch()
% This function should be deleted later. It's just for holding the
% reference code

% I'll add a list of some of the polynomials to use for 
% checking results


% Hermite 

h0 = 1; h1 = [1 0]; h2 = [1 0 -1]; h3 = [1 0 -3 0]; 
h4 = [1 0 -6 0 3]; h5 = [1 0 -10 0 15 0];
h6 = [1 0 -15 0 45 0 -15]; h7 = [1 0 -21 0 105 0 -105];
h8 = [1 0 -28 0 210 0 -420 0 105];
h9 = [1 0 -36 0 378 0 1260 0 945 0];
h10 = [1 0 -45 0 630 0 -3150 0 4725 0 -945];

% Legendre
l0 = 1; l1 = [1 0]; l2 = [3 0 -1]/2; l3 = [5 0 -3 0]/2;
l4 = [35 0 -30 0 3]/8; l5 = [63 0 -70 0 15 0]/8;
l6 = [231 0 -315 0 105 -5]/16; l7 = [429 0 -693 0 315 0 -35 0]/16;
l8 = [6435 0 -12012 0 6930 0 -1260 0 35]/128;
l9 = [12155 0 -25740 0 18018 0 -4620 0 315 0]/128;
l10 = [46189 0 -109395 0 90090 0 -30030 0 3465 0 -63]/256;

% Order to use
N = 4; 


% Hermite recurrence relation (Gaussian weight)
a = ones(N,1);
b = zeros(N,1);
c = (0:(N-1))';

% Legendre recurrence relation (uniform weight)
n = N-1; 
a = (2*(0:n) +1)'./((0:n)+1)';
b = zeros(N,1);
c = (0:n)'./((0:n)+1)';

% Jacobi recurrence relation (beta weight)
% These are the alpha and beta parameters for the beta distribution
al = 3; bt = 2;
al = 0; bt = 0; 
% Denominator of this ridiculous expression
n = (0:(N-1)); 
a = 0.5 * (2*n+al+bt+1) .* (2*n+al+bt+2) ...
        ./ ((n+1) .* (n+al+bt+1));
b = 0.5 * (al^2-bt^2) .* (2*n+al+bt+1) ...
        ./ ((n+1) .* (n+al+bt+1) .* (2*n+al+bt));
c =       (n+al).*(n+bt).*(2*n+al+bt+2) ...
        ./ ((n+1) .* (n+al+bt+1) .* (2*n+al+bt));

% Correction to first term
a(1) = 0.5 * (al + bt + 2);
b(1) = 0.5 * (al - bt); 


 

end
