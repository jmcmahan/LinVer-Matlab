function sample = draw_posterior_sample(param, post, M)
% draw_posterior_sample(param, post, N) - Draw N samples from the posterior
% defined by post for the problem defined by param.

if strcmp(param.unknowns, 'beta')
    L = chol(post.sigma2 / param.lambda);
    sample = zeros(M, param.Nbeta);
    for j = 1:M
        sample(j, :) = post.mu2' + randn(1, param.Nbeta)*L;
    end
elseif strcmp(param.unknowns, 'beta_lambda')
    L = chol(post.sigma2);
    sample = zeros(M, param.Nbeta+1);
    for j = 1:M
        l = gamrnd(post.a1, 1/post.b1);
        sample(j, 1:param.Nbeta) = post.loc' + randn(1, param.Nbeta)*L/sqrt(l);
        sample(j, end) = l;
    end
elseif strcmp(param.unknowns, 'beta_lambda_phi')
    sample = zeros(M, param.Nbeta+2);
    phi = phirnd(post, param, M);
    % Reuse the code for evaluating the lambda / beta parameters with the
    % samples of phi drawn here.
    paramj = param;
    paramj.unknowns = 'beta_lambda';
    for j = 1:M
        paramj.phi = phi(j);
        postj = eval_posterior(paramj);
        L = chol(postj.sigma2);
        l = gamrnd(postj.a1, 1/postj.b1);
        sample(j, 1:paramj.Nbeta) = postj.loc' + randn(1, paramj.Nbeta)*L/sqrt(l);
        sample(j, end-1) = l;
        sample(j, end) = phi(j);
    end
end

end
