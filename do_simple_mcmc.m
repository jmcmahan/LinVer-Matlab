function [qchain, aratio, oob] = do_simple_mcmc(param, post, Nchain)
% [chain, aratio] = do_simple_mcmc(param) - Use a Metropolis-Hastings algorithm to generate
% a 10,000 iterate chain drawn from the posterior of the problem defined by param. 
% Acceptance ratio returned as aratio.
% 
% [chain, aratio] = do_simple_mcmc(param, post) - Same as above but use the posterior 
% estimates to initialize the proposal matrix.
%
% [chain, aratio] = do_simple_mcmc(param, post, Nchain) - Same as above, but generate Nchain 
% samples.
%
% NOTE: This does not use any of the unknown parameters in anyway. 

if nargin < 3
    Nchain = 1e4; 
end

Nbeta = param.Nbeta;
y = param.y;
G = param.G;

if strcmp(param.unknowns, 'beta')
    calcase = 1;
    Np = Nbeta;
    ll = @(q) log_likelihood(param, q);
elseif strcmp(param.unknowns, 'beta_lambda')
    calcase = 2;
    Np = Nbeta + 1;
    ll = @(q) log_likelihood(param, q(1:Nbeta), q(Nbeta+1));
elseif strcmp(param.unknowns, 'beta_lambda_phi')
    calcase = 3;
    Np = Nbeta + 2;
    ll = @(q) log_likelihood(param, q(1:Nbeta), q(Nbeta+1), q(Nbeta+2));
end


% Set up parameter ranges
if calcase == 1
    qrange = param.betarange;
elseif calcase == 2
    qrange = [param.betarange; param.lambdarange];
elseif calcase == 3
    qrange = [param.betarange; param.lambdarange; param.phirange];
end


qchain = zeros(Nchain, Np);
% For now, initialize to the true values. This isn't as unfair as it 
% sounds since we're looking at distributions, not the single guess, 
% and this wouldn't guarantee the distribution is correct. 

qchain(1, 1:Nbeta) = param.beta;
if calcase > 1
    qchain(1, Nbeta+1) = param.lambda;
end
if calcase > 3
    qchain(1, Nbeta+2) = param.phi;
end

%res = y - G*post.mu2;
res = y - G*param.beta;
s2ols = res'*res / (param.N - Np);
Ri = eval_corrfuncinv(param);
%V = s2ols*param.lambda*inv(G'*Ri*G);
V = zeros(Np);
V(1:param.Nbeta, 1:param.Nbeta) = s2ols*inv(G'*Ri*G);
if calcase > 1
    % For now, just put something reasonable for the range of the hyper
    % parameters we're using in testing. 
    %V(param.Nbeta+1, param.Nbeta+1) = 10;
    lrange = param.lambdarange(2) - param.lambdarange(1);
    V(param.Nbeta+1, param.Nbeta+1) = lrange / 100;
end
if calcase > 2
    prange = param.phirange(2) - param.phirange(1);
    V(end, end) = prange / 100;
end

L = chol(V)';

acceptnum = 1;
oob = 0;


for k = 2:Nchain
    accept = false;
    if ~mod(k, fix(Nchain/100))
        disp(sprintf('Percent complete: %d%%\tAcceptance Ratio: %d', ...
                        round(100*k/Nchain), acceptnum/(k-1)))
    end

    qp = qchain(k-1, :)';
    qstar = qp + L*randn(Np, 1);

    if sum(qstar <= qrange(:, 1)) || sum(qstar >= qrange(:, 2))
        % Set the probability to 0 if the sample is out-of-bounds
        oob = oob + 1; 
        r = 0; 
    else
        r = exp(ll(qstar) - ll(qp));
    end

    if r >= 1
        qchain(k,:) = qstar';
        accept = true;
    else
        flip = rand;
        if flip < r
            qchain(k,:) = qstar';
            accept = true;
        else
            qchain(k,:) = qp';
        end
    end
    if accept
        acceptnum  = acceptnum + 1;
    end
end


aratio = acceptnum/Nchain;

end





function ll = log_likelihood(param, beta, lambda, phi)
% Evaluate the log-likelihood, using the number of arguments passed to 
% determine the correct form of the likelihood function to use.

G = param.G; y = param.y;
N = param.N; Nbeta = param.Nbeta;
prior = param.prior;

if nargin < 4
    phi = param.phi;
end
if nargin < 3
    lambda = param.lambda;
end

Ri = eval_corrfuncinv(param, phi);

res = y - G*beta;

if strcmp(prior.type, 'noninformative')
    if nargin == 2
        % Beta unknown
        ll = - 0.5 * lambda * res'*Ri*res;  
    elseif nargin == 3
        % Beta, lambda unknown
        ll = (0.5*N-1) * log(lambda) ...    
            - 0.5 * lambda * res'*Ri*res;   
    elseif nargin == 4
        % Beta, lambda, phi unknown
        d = eval_det(param, phi); 
        ll = -  0.5*log(d) ...              
             + (0.5*N-1) * log(lambda) ...  
             -  0.5 * lambda * res'*Ri*res; 
    end
elseif strcmp(prior.type, 'gaussian')
    mu0 = prior.mu0;    sigma0 = prior.sigma0;
    resbeta = mu0 - beta;
    if nargin == 2
        % Beta unknown
        ll = - 0.5 * lambda * res'*Ri*res   
             - 0.5 * lambda * resbeta'*inv(sigma0)*resbeta;
    elseif nargin == 3
        % Beta, lambda unknown
        ll = (0.5*(N+Nbeta) - 1) * log(lambda) ...
            - 0.5 * lambda * res'*Ri*res
            - 0.5 * lambda * resbeta'*inv(sigma0)*resbeta;
    elseif nargin == 4
        % Beta, lambda, phi unknown
        d = eval_det(param, phi); 
        ll = -  0.5*log(d) ...
             + (0.5*(N+Nbeta)-1) * log(lambda) ...
             -  0.5 * res'*Ri*res / lambda;
             - 0.5 * lambda * resbeta'*inv(sigma0)*resbeta;
    end
end

end
