function [qchain, aratio, oob] = do_simple_mcmc_lambda(param, post, Nchain)
% [chain, aratio] = do_simple_mcmc_lambda(param) - Use a Metropolis-Hastings algorithm to generate
% a 10,000 iterate chain drawn from the posterior of the problem defined by param. 
% Acceptance ratio returned as aratio.
% 
% [chain, aratio] = do_simple_mcmc_lambda(param, post) - Same as above but use the posterior 
% estimates to initialize the proposal matrix.
%
% [chain, aratio] = do_simple_mcmc_lambda(param, post, Nchain) - Same as above, but generate Nchain 
% samples.
%
% This code is similar to do_simple_mcmc, but samples lambda from a gamma
% distribution rather than including it as part of the Gaussian proposal
% matrix. 
%
%
% NOTE: This does not use any of the unknown parameters in anyway. 

if nargin < 3
    Nchain = 1e4; 
end

Nbeta = param.Nbeta;
y = param.y;
G = param.G;
N = param.N;

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


% NOTE: Using the actual value of beta for simplicity. A somewhat better
% exercise would be to use the MLE estimate. Since we are verifying the
% distribution, though, and not just the point estimate, this should still
% provide a reasonable test. 
% Also worth noting - this is using the true value of phi to evaluate
% Ri. This is also additional accuracy, but again, shouldn't be too
% detrimental for this test. 

res = y - G*param.beta;
Ri = eval_corrfuncinv(param);
param.Ri = Ri;
s2ols = res'*Ri*res / (param.N - Np);

% Lambda will be sampled from a gamma distribution rather than the
% Gaussian proposal.
if calcase == 1
    V = zeros(Np);
else
    V = zeros(Np-1);
end

if calcase > 1
    V(1:param.Nbeta, 1:param.Nbeta) = inv(G'*Ri*G) * s2ols;
else
    % Case 1 has lambda known, so use it here.
    V(1:param.Nbeta, 1:param.Nbeta) = inv(G'*Ri*G) / param.lambda;
end


if calcase > 1
    % For now, just put something reasonable for the range of the hyper
    % parameters we're using in testing. 
    %V(param.Nbeta+1, param.Nbeta+1) = 10;
    lrange = param.lambdarange(2) - param.lambdarange(1);
end

if calcase > 2
    prange = param.phirange(2) - param.phirange(1);
    V(end, end) = prange / 100;
end

L = chol(V)';

acceptnum = 1;
oob = 0;

% Parameters for sampling lambda from the gamma distribution 
Ns = 0.01;
s2s = s2ols; 
alpha = 0.5 * (Ns + N); 

for k = 2:Nchain
    accept = false;
    if ~mod(k, fix(Nchain/100))
        disp(sprintf('Percent complete: %d%%\tAcceptance Ratio: %d', ...
                        round(100*k/Nchain), acceptnum/(k-1)))
    end

    qp = qchain(k-1, :)';
    if calcase == 1
        qstar = qp + L*randn(length(L), 1);
    elseif calcase == 2
        qstar = qp(1:Nbeta) + L*randn(length(L), 1);
        % Add in lambda, sampled from the gamma distribution
        qstar = [qstar(1:Nbeta); qp(end)]; 
    elseif calcase == 3
        qstar = qp([1:Nbeta, Nbeta+2]) + L*randn(length(L), 1);
        % Add in in lambda, sampled from gamma distribution
        qstar = [qstar(1:Nbeta); qp(end-1); qstar(end)]; 
    end
    
    if calcase > 1
        if calcase == 2
            ss = sum_of_squares(param, qstar(1:end-1));
        elseif calcase == 3
            ss = sum_of_squares(param, qstar(1:end-2), qstar(end));
        else
            disp('ERROR: Invalid calibration cases (should never see this).')
        end
        
        beta = 0.5 * (Ns * s2s + ss);
        % New lambda sample. Note this is put in the chain afterwards,
        % regardless of rejection / acceptance of the proposed parameters.
        lnew = gamrnd(alpha, 1 / beta); 
    end
    
    
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
    
    if calcase == 2
        qchain(k, end) = lnew;
    elseif calcase == 3
        qchain(k, end-1) = lnew;
    end
end


aratio = acceptnum/Nchain;

end


function ss = sum_of_squares(param, beta, phi)
    if nargin < 3
%        Ri = eval_corrfuncinv(param);
        % Abuse! We've added Ri in this case. When phi is known this
        % saves from doing a Cholesky factorization each iteration
        % of the chain. 
        Ri = param.Ri;
    else
        Ri = eval_corrfuncinv(param, phi);
    end
    if isrow(beta)
        beta = beta';
    end
    G = param.G;
    y = param.y;
    res = y - G*beta; 
    ss = res' * Ri * res;
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
