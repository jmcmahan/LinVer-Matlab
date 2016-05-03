function qchain = do_simple_mcmc(param, Nchain)
% chain = do_simple_mcmc(param) - Use a Metropolis-Hastings algorithm to generate
% a 10,000 iterate chain drawn from the posterior of the problem defined by param.
% 
% chain = do_simple_mcmc(param, Nchain) - Same as above, but generate Nchain 
% samples.
%
% NOTE: This does not use any of the unknown parameters in anyway. 

if nargin < 2
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

res = y - G*param.beta;
s2ols = res'*res / (param.N - Np);
Ri = eval_corrfuncinv(param);
V = s2ols*param.lambda*inv(G'*Ri*G);
L = chol(V)';


accrej = zeros(Nchain, 1);
for k = 2:Nchain
    if ~mod(k, fix(Nchain/100))
        disp(sprintf('Percent complete: %d%%', round(100*k/Nchain)))
    end

    qp = qchain(k-1, :)';
    qstar = qp + L*randn(Np, 1);
    % Add something to ensure samples are in-bounds here

    r = exp(ll(qstar) - ll(qp));
    if r >= 1
        qchain(k,:) = qstar';
        accrej(k) = 1; 
    else
        flip = rand;
        if flip < r
            qchain(k,:) = qstar';
            accrej(k) = 1;
        else
            qchain(k,:) = qp';
        end
    end
end



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
        d = eval_det(phi, param); 
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
        d = eval_det(phi, param); 
        ll = -  0.5*log(d) ...
             + (0.5*(N+Nbeta)-1) * log(lambda) ...
             -  0.5 * res'*Ri*res / lambda;
             - 0.5 * lambda * resbeta'*inv(sigma0)*resbeta;
    end
end

end
