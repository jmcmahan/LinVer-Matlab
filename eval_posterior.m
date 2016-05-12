function post = eval_posterior(param)
% post = eval_posterior(param) - Evaluate the posterior parameters for the 
% problem specified by the structure "param". 

global debug;

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

post.debug = debug;

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
    [post.pbeta, post.plambda, post.pphi, post.xphi, post.c, ...
        post.phiparam]  ...
          = eval_beta_lambda_phi_params(param);
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


function [a1, b1, dof, loc, scl, sigma2, sigma3, GRiG] = eval_beta_lambda_params(param, phi)
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
% GRiG = another term used when phi is unknown

    prior = param.prior;
    G = param.G;    y = param.y;
    Ri = eval_corrfuncinv(param, phi); 
    N = param.N;    Nbeta = param.Nbeta;

    GRiG = G'*Ri*G; 
    if strcmp(prior.type, 'noninformative')
        % Non-informative prior
        sigma2 = inv(GRiG);
        mu2 = sigma2*(G'*Ri*y);
        betahat = mu2; 
        % Gamma parameters for lambda 
        a1 = (N-Nbeta)/2;
        resy = y - G*betahat; 
        b1 = (resy'*Ri*resy) / 2;
        sigma3 = [];
    elseif strcmp(prior.type, 'gaussian')
        % Gaussian prior
        sigma0 = prior.sigma0; mu0 = prior.mu0;
        %sigma2 = inv(inv(sigma0) + G'*Ri*G);
        sigma2 = inv(inv(sigma0) + GRiG);
        mu2 = sigma2 * (G'*Ri*y + sigma0\mu0);
        betahat = (GRiG) \ (G'*Ri*y);
        sigma3 = sigma0 + inv(GRiG); 
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


function [pbeta, plambda, pphi, xphi, c, phiparam] = eval_beta_lambda_phi_params(param)
% Find the posterior distribution parameters when calibrating beta, lambda, and
% phi. 
%
% pbeta - Function to evaluate marginal beta posterior at beta. So do pbeta(beta)
%         to get the probability at beta. 
%
% plambda - Same as above but for lambda.
%
% pphi, xphi - Unlike the above, this is not a function, but a set of nodes
%              and the values at those nodes. 
%
% c - Integration constant for normalizing the phi posterior. 
% 
% phiparam - structure containing evaluation of the parameters for the lambda and 
%            phi posteriors at the points xphi

    global debug;

    prior = param.prior;
    phirange = param.phirange;
    lambdarange = param.lambdarange;
    betarange = param.betarange;
    Nbeta = param.Nbeta;
    phitrue = param.phi;

    % This determines the prior for phi
    if ~isfield(prior, 'phitype')
        phitype = 'uniform';
    else
        phitype = prior.phitype;
    end

    % Order of the quadrature used for computing phi
    if ~isfield(param, 'quadorder')
        quadorder = 100;
    else
        quadorder = param.quadorder;
    end

    if strcmp(param.corrfunc, 'none')
        disp('Error in posterior evaluation: Cannot estimate phi for the uncorrelated correlation function');
        return
    end

    % Nodes and weights for the appropriate quadrature
    [xphi, wphi] = gauss_quadrature(phitype, quadorder);

    % Map the nodes onto the domain of interest
    if strcmp(phitype, 'uniform') || strcmp(phitype, 'beta')
        a = phirange(1); b = phirange(2);
        xphi = (xphi + 1)*(b - a)/2 + a; 
    elseif strcmp(phitype, 'gaussian')
        % Not currently implemented - this needs a mean and variance to be defined
        % for the phi prior 
        %xphi = xphi*sigma + mu;
    end


    disp('Computing phi...')
    measure.alpha = param.N;
    measure.beta = param.N;
    [pphi, a1, b1, dof, scl, loc] = phi_integrand(param, xphi, measure);
    c = wphi'*pphi;
    pphi = pphi / c; 

    % This is all used to build the functions for evaluating the marginal beta and lambda
    phiparam.quadorder = quadorder;
    phiparam.pphi = pphi; 
    phiparam.xphi = xphi;
    phiparam.wphi = wphi;
    phiparam.Nbeta = Nbeta;
    phiparam.a1 = a1; phiparam.b1 = b1;
    phiparam.dof = dof; phiparam.scl = scl; phiparam.loc = loc; 

    if strcmp(phitype, 'uniform')
        % The prior term of the phi posterior is separated into the quadrature weights (wphi)
        % so for plotting purposes, it needs to be incorporated back in. Only have this
        % for the uniform prior since that seems to be working fine, but this needs to be
        % done for the other phi priors if they are implemented
        pphi = pphi / (b-a);
    end
    disp('Done.');
    

    % This uses quadrature nodes to determine which values of beta and lambda to 
    % evaluate.

    pbeta = @(beta) beta_marginal(phiparam, beta); 
    plambda = @(lambda) lambda_marginal(phiparam, lambda);


    debug.wphi = wphi;
end


function [f, a1, b1, dof, scl, loc] = phi_integrand(param, phi)
% f = phi_integrand(param, post, phi) : Evaluate part of the integrand involved
% in the calculation of the posterior of phi. The problem is defined
% in param and phi is a vector of values at which to evaluate the
% integrand. Note that the values computed here do not include the value of 
% the prior of phi. 


G = param.G;
M = length(phi);
Nbeta = param.Nbeta;

f = zeros(size(phi));
a1 = zeros(size(phi));
b1 = zeros(size(phi));
dof = zeros(size(phi));
scl = zeros(length(phi), Nbeta, Nbeta);
loc = zeros(length(phi), Nbeta);

% This computation can be done outside of the loop, for whatever good that does
detR = eval_det(param, phi);
for j = 1:M
    phic = phi(j); 
    [a1j, b1j, dofj, locj, sclj, ~, sigma3, GRiG] = eval_beta_lambda_params(param, phic);

    detRc = detR(j); 
    % NOTE - the code as-is returns sigma3 = [] in the non-informative case. 
    % In Matlab, det([]) == 1, so this should be fine.
    dets3 = det(sigma3);
    detGRiG = det(GRiG);

    f(j) = 1 / (b1j^a1j * sqrt(detRc) * sqrt(detGRiG) * sqrt(dets3));

    a1(j) = a1j; b1(j) = b1j; 
    dof(j) = dofj; scl(j,:,:) = sclj; loc(j,:) = locj';
end

end


function p = beta_marginal(phiparam, beta)
% This is the marginal posterior for beta when phi is unknown. This evaluates
% the probability at beta given info in phiparam. The parameter beta should be
% size M by Nbeta where Nbeta is the number of regression parameters (dimension
% of the linear problem) and M is the number of points to evaluate p(beta) at.


    quadorder = phiparam.quadorder;
    pphi = phiparam.pphi;
    xphi = phiparam.xphi;
    wphi = phiparam.wphi;
    Nbeta = phiparam.Nbeta;
    a1 = phiparam.a1; b1 = phiparam.b1;
    dof = phiparam.dof; scl = phiparam.scl; loc = phiparam.loc; 


    if Nbeta == 1
        Ns = length(beta);
        if ~isrow(beta)
            beta = beta';
        end

        ignd = zeros(quadorder, Ns);
        % Can vectorize this, as well, if a speed-up is needed at some point
        for j = 1:quadorder
            ignd(j,:) = tpdf( (beta - loc(j)) / sqrt(scl(j)), dof(j)) ...
                                / sqrt(scl(j));
        end
        p = (wphi.*pphi)'*ignd;
    else
        disp('Case 3 currently only supports 1 regression parameter');
        return;
    end
    
end


function [p, ignd] = lambda_marginal(phiparam, lambda)
% This is the marginal posterior for lambda when phi is unknown. This evaluates
% the probability at lambda given info in phiparam. The parameter beta should be
% size M by 1 where M is the number of points to evaluate p(lambda) at.


    quadorder = phiparam.quadorder;
    pphi = phiparam.pphi;
    xphi = phiparam.xphi;
    wphi = phiparam.wphi;
    Nbeta = phiparam.Nbeta;
    a1 = phiparam.a1; b1 = phiparam.b1;
    dof = phiparam.dof; scl = phiparam.scl; loc = phiparam.loc; 


    Ns = length(lambda);
    if ~isrow(lambda)
        lambda = lambda';
    end

    ignd = zeros(quadorder, Ns); 
    for j = 1:quadorder;
        ignd(j, :) = gampdf(lambda, a1(j), 1 / b1(j));
    end

    p = (wphi.*pphi)'*ignd;
    
end



