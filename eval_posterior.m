function post = eval_posterior(param)
% post = eval_posterior(param) - Evaluate the posterior parameters for the 
% problem specified by the structure "param". 


corrfunc = param.corrfunc;
unknowns = param.unknowns;


if strcmp(unknowns, 'beta')
    post = do_beta(param);
elseif strcmp(unknowns, 'beta_lambda')
    post = do_beta_lambda(param);
elseif strcmp(unknowns, 'beta_lambda_phi')
    post = do_beta_lambda_phi(param);
else
    disp('Error in eval_posterior: Unsupported unknowns')
end


end



function post = do_beta(param)
    phi = param.phi;
    [post.sigma2, post.mu2] = eval_beta_params(param, phi);
end


function post = do_beta_lambda(param)
    phi = param.phi;
    [post.a1, post.b1, post.dof, post.loc, post.scl, post.sigma2] = ...
                eval_beta_lambda_params(param, phi); 
end


function post = do_beta_lambda_phi(param)
end


function [sigma2, mu2] = eval_beta_params(param, phi)
    % Find the posterior distribution parameters when calibrating beta.
    %
    % sigma2 = unscaled posterior covariance matrix of beta
    % mu2 = posterior mean of beta
    prior = param.prior;
    G = param.G;    y = param.y;
    Ri = eval_corrfuncinv(param, phi);

    if strcmp(prior.type, 'noninformative')
        % Non-informative prior
        sigma2 = inv(G'*Ri*G);
        mu2 = sigma2*(G'*Ri*y);
    elseif strcmp(prior.type, 'gaussian')
        sigma0 = prior.sigma0;  mu0 = prior.mu0;
        sigma2 = inv(inv(sigma0) + G'*Ri*G);
        mu2 = sigma2 * (G'*Ri*y + sigma0\mu0);
    end
end


function [a1, b1, dof, loc, scl, sigma2, sigma3] = eval_beta_lambda_params(param, phi)
    % Find the posterior distribution parameters when calibrating beta and 
    % lambda.
    % 
    % a1, b1 = Gamma distribution parameters. 
    % NOTE: The notation of parameters for the Gamma distribution varies. 
    % Matlab uses a = a_1, b = 1/b1 for its "gampdf" function. 
    %
    % dof = degrees of freedom for t-distribution of beta posterior
    % loc = location parameter for t-distribution of beta posterior
    % scl = scale parameter for t-distribution of beta posterior
    % NOTE: In Matlab, the density for a 1D t-distribution is 
    % (1/sqrt(scl)) * tpdf((t-loc)/sqrt(scl), dof)
    % 
    % sigma3 = additional covariance matrix used when phi is unknown

    prior = param.prior;
    G = param.G;    y = param.y;
    Ri = eval_corrfuncinv(param, phi); 
    N = param.N;    Nbeta = param.Nbeta;

    if strcmp(prior.type, 'noninformative')
        % Non-informative prior
        sigma2 = inv(G'*Ri*G);
        mu2 = sigma2*(G'*Ri*y);
        betahat = mu2; 
        % Gamma parameters for lambda 
        a1 = (N-Nbeta)/2;
        resy = y - G*betahat; 
        b1 = (resy'*Ri*resy) / 2;
    elseif strcmp(prior.type, 'gaussian')
        % Gaussian prior
        sigma0 = prior.sigma0; mu0 = prior.mu0;
        sigma2 = inv(inv(sigma0) + G'*Ri*G);
        mu2 = sigma2 * (G'*Ri*y + sigma0\mu0);
        betahat = (G'*Ri*G) \ (G'*Ri*y);
        sigma3 = sigma0 + inv(G'*Ri*G); 
        % Gamma parameters for lambda
        a1 = N / 2;
        resy = y - G*betahat;
        resbeta = betahat - mu0;
        b1 = (resy'*Ri*resy + resbeta'*(sigma3 \ resbeta));
    end

    

    % t-distribution parameters for beta
    dof = 2*a1;
    loc = mu2;
    scl = b1*sigma2/a1;

end
