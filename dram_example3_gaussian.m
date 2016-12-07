% Equicorrelation, beta, lambda fit
    param.N = 300;                          % Number of data points
    param.beta = [1.5; 3.5];                 % True regression parameters
    param.Nbeta = length(param.beta);       % Number of regression parameters
    param.lambda = 10.0;                    % True scale parameter 
    param.phi = 0.2;                        % True correlation parameter
    param.corrfunc = 'equal';               % Correlation function 
    param.prior.type = 'gaussian';    % Use non-informative prior
    param.unknowns = 'beta_lambda';         % Compute posterior for the regression parameters only
    param.G = randn(param.N, param.Nbeta);  % Design matrix...
    param.G(:,1) = 1;                       % ...with first regression parameter a bias term.
    param.betarange = [-100*ones(param.Nbeta, 1), 100*ones(param.Nbeta, 1)];
    param.lambdarange = [1e-1, 1e3];
    param.prior.sigma0 = 0.1*eye(param.Nbeta);        % Prior unscaled covariance
    param.prior.mu0 = [2; 3];             % Prior mean

    y0 = param.G * param.beta;              % Error-free observation data

    e = eval_noise(param);                  % Generate unscaled observation error
    param.y = y0 + e / sqrt(param.lambda);  % Add error to create data for calibration
    y = param.y;
    G = param.G;
    post = eval_posterior(param);           % Compute posterior

    %br = linspace(-1, 5)';                  % Range to evaluate the beta posterior on
    b1 = linspace(post.loc(1) - sqrt(post.scl(1,1))*4.5, post.loc(1) + sqrt(post.scl(1,1))*4.5, 500);
    b2 = linspace(post.loc(2) - sqrt(post.scl(2,2))*4.5, post.loc(2) + sqrt(post.scl(2,2))*4.5, 500);
    [brx, bry] = meshgrid(b1,b2);
    br = [brx(:), bry(:)];
    bpost = post.pbeta(br);

    l = linspace(param.lambdarange(1), param.lambdarange(end), 1e4)';
    lmu = post.a1 / post.b1;
    lsd = sqrt(post.a1 / post.b1^2);
    l = linspace(max(lmu - 4*lsd, 0.1) , lmu + 4*lsd, 1e3);
    lpost = post.plambda(l); 

    %plot(br, bpost);               % Plot beta posterior

 
    % 100 samples from a distribution different than the posterior
    bad_sample = [randn(100,2), exp(randn(100,1))];

    % 100 samples from same distribution as posterior

    % This would be a less manual way to draw from the posterior
    good_sample = draw_posterior_sample(param, post, 100); 

    param.nsimu = 100000;

    % Draw from DRAM to verify
    sol = dram_from_linver(param, false);

    good_sample = sol.chain(20001:500:end, :);

    % Draw sample from DRAM with bug added to sum-of-squares function

    solbad = dram_from_linver(param, true);
    bad_sample = solbad.chain(20001:500:end,:);


    % Confidence level
    alpha = 0.01;
    % Number of tests to run
    numtests = 500;

    % Do the energy tests for each of these samples
    disp('Test Key:')
    disp('1 Indicates Sample Significantly Different.')
    disp('0 Indicates Sample No Significant Difference Found.')
    result_good = do_energy_test(good_sample, param, post, alpha, numtests);
    result_bad = do_energy_test(bad_sample, param, post, alpha, numtests); 

    disp(' ')
    disp('Good sample results:')
    disp(result_good.fail_ratio) 
    disp('Should not significantly exceed alpha = 0.1.')
    disp(' ')
    disp('Bad sample results:')
    disp(result_bad.fail_ratio)
    disp('Hopefully exceeds alpha = 0.1. ')
    disp(' ')
    disp('Note that small sample size used may allow number of failures')
    disp('in the good sample to be somewhat higher than 0.1. Ratio should')
    disp('be much smaller than that of the bad sample, though.')

    save('dram_example3_gaussian.mat', 'result_good', 'result_bad', 'good_sample', 'bad_sample', 'sol', 'solbad', 'G', 'y', 'e', 'post', '-mat')
