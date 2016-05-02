function chain = do_simple_mcmc(param, Nchain)
% chain = do_simple_mcmc(param) - Use a Metropolis-Hastings algorithm to generate
% a 10,000 iterate chain drawn from the posterior of the problem defined by param.
% 
% chain = do_simple_mcmc(param, Nchain) - Same as above, but generate Nchain 
% samples.
%
% NOTE: This does not use any of the unknown parameters in anyway. 




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
