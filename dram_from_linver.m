function sol = dram_from_linver(param)
    G = param.G;
    Nbeta = param.Nbeta;
    
    % We'll set this one up to do all the parameter inference in the 
    % chain, i.e., without using the gamma distribution to sample 
    % lambda
    
    
    Cinv = eval_corrfuncinv(param);
    
    if strcmp(param.unknowns, 'beta')
        % Lambda is known in this case, so use it
        qcov = inv(G'*Cinv*G) / param.lambda; 
        Nparam = Nbeta; 
        
    elseif strcmp(param.unknowns, 'beta_lambda')
        % This is using the actual parameters for initialization,
        % since it's simpler. It should be ok, as we're verifying the
        % distribution, not any particular point estimate. The idea is
        % that we are assuming this can be initialized in some reasonable
        % way that gives similar results. 
        res = param.y - G*param.beta; 
        s2ols = res'*Cinv*res / (param.N - Nbeta - 1); 
        qcov = zeros(Nbeta+1);
        % Beta part of the initial proposal covariance
        qcov(1:Nbeta,1:Nbeta) = inv(G'*Cinv*G) * s2ols;
        % Lambda part of the initial proposal covariance
        qcov(Nbeta+1,Nbeta+1) = 10;
        Nparam = Nbeta + 1;
        
    elseif strcmp(param.unknowns, 'beta_lambda_phi')
        res = param.y - G*param.beta; 
        % Note this is also using the true information about phi, since
        % Cinv is computed with the actual parameter. 
        s2ols = res'*Cinv*res / (param.N - Nbeta - 1); 
        qcov = zeros(Nbeta+2);
        % Beta part of the initial proposal covariance
        qcov(1:Nbeta,1:Nbeta) = inv(G'*Cinv*G) * s2ols;
        % Lambda part of the initial proposal covariance
        qcov(Nbeta+1,Nbeta+1) = 10;  
        % Phi part of the initial proposal covariance
        qcov(Nbeta+2,Nbeta+2) = 5e-2;
        Nparam = Nbeta + 2; 
        
    end
    

    
    drampar = cell(Nparam, 1);
    for j = 1:Nbeta
        if strcmp(param.prior.type, 'noninformative')
            drampar{j} = {'', param.beta(j), param.betarange(j, 1), ...
                        param.betarange(j, 2)};
        elseif strcmp(param.prior.type, 'gaussian')
            drampar{j} = {'', param.beta(j), param.betarange(j, 1), ...
                        param.betarange(j, 2), ...
                        param.prior.mu0(j), ...
                        sqrt(param.prior.sigma0(j,j))};            
        end
    end
    if strcmp(param.unknowns, 'beta_lambda') ...
            || strcmp(param.unknowns, 'beta_lambda_phi')
        drampar{Nbeta+1} = {'', param.lambda, param.lambdarange(1), ...
                        param.lambdarange(2)}; 
    end
    if strcmp(param.unknowns, 'beta_lambda_phi')
        drampar{end} = {'', param.phi, param.phirange(1), ...
                            param.phirange(2)};
    end
 
    if strcmp(param.unknowns, 'beta')
        % Use the known true-value of lambda
        priorfun = @(th, mu, sig) param.lambda * sum( ((th-mu)./sig).^2 );
        
    else
        % True value of lambda unknown, so use the sampled value
        priorfun = @(th, mu, sig) -param.Nbeta*log(th(Nbeta+1)) + th(Nbeta+1) * ...
                     sum( ((th(1:Nbeta)-mu(1:Nbeta)) ./sig(1:Nbeta)).^2 );
                 
    end

    data.ydata = param.y';
    data.xdata = linspace(0, 1, param.N);  
    
    % The model sigma is set to 1, since we are trying to infer it as
    % part of the random walk. May be interesting to try this with sigma
    % updated automatically.
    model.sigma2 = 1; 
    model.priorfun = priorfun;
    model.N = param.N;
    if strcmp(param.unknowns, 'beta')
        model.ssfun = @(theta, data) ssfun1(theta, data, param, Cinv);
    elseif strcmp(param.unknowns, 'beta_lambda')
        model.ssfun = @(theta, data) ssfun2(theta, data, param, Cinv);
    elseif strcmp(param.unknowns, 'beta_lambda_phi')
        model.ssfun = @(theta, data) ssfun3(theta, data, param);
    end
    options.qcov = qcov;
    options.updatesigma = 0;
    
    [res, chain] = mcmcrun(model, data, drampar, options);
    sol.res = res;
    sol.chain = chain;
    
end



function ss = ssfun1(theta, data, param, Cinv)
    % Log-likelihood function when beta is unknown
    G = param.G;
    if isrow(theta); 
        theta = theta';
    end
    
    if isrow(data.ydata)
        y = data.ydata';
    end
    
    yhat = G*theta;
    resid = y - yhat;
    ss = resid'*Cinv*resid * param.lambda; 
end

function ss = ssfun2(theta, data, param, Cinv)
    % Log-likelihood function when beta, lambda are unknown
    G = param.G;
    if isrow(theta); 
        theta = theta';
    end
    
    if isrow(data.ydata)
        y = data.ydata';
    end
    
    beta = theta(1:param.Nbeta);
    lambda = theta(param.Nbeta+1);
    yhat = G*beta;
    resid = y - yhat;
    ss = resid'*Cinv*resid * lambda - param.N * log(lambda);
end

function ss = ssfun3(theta, data, param)
    % Log-likelihood function when beta,lambda,phi are unknown    
    G = param.G;
    if isrow(theta); 
        theta = theta';
    end
    
    if isrow(data.ydata)
        y = data.ydata';
    end
    
    beta = theta(1:param.Nbeta);
    lambda = theta(param.Nbeta+1);
    phi = theta(param.Nbeta+2);
    Cinv = eval_corrfuncinv(param, phi); 
    detC = eval_det(param, phi);
    
    yhat = G*beta;
    resid = y - yhat;
    
    ss = resid'*Cinv*resid * lambda + log(detC) - param.N * log(lambda);
end

