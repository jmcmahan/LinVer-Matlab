% AR(1) correlation, beta, lambda, phi fit
    %param.N = 300;                          % Number of data points
    param.N = 300;                          % Number of data points
    param.beta = [1.5; 3.5];                 % True regression parameters
    param.Nbeta = length(param.beta);       % Number of regression parameters
    param.lambda = 10.0;                    % True scale parameter 
    param.phi = 0.5;                        % True correlation parameter
    param.corrfunc = 'ar';               % Correlation function 
    param.prior.type = 'noninformative';    % Use non-informative prior
    param.unknowns = 'beta_lambda_phi';     % Compute posterior for the regression parameters only
    param.G = randn(param.N, param.Nbeta);  % Design matrix...
    param.G(:,1) = 1;                       % ...with first regression parameter a bias term.
    param.betarange = [-100*ones(param.Nbeta, 1), 100*ones(param.Nbeta, 1)];
    param.lambdarange = [1e-1, 1e3];
    param.phirange = [-0.95, 0.95];
    
    y0 = param.G * param.beta;              % Error-free observation data

    e = eval_noise(param);                  % Generate unscaled observation error
    param.y = y0 + e / sqrt(param.lambda);  % Add error to create data for calibration
    y = param.y;
    G = param.G;
    post = eval_posterior(param);           % Compute posterior

    if 1 == 0 
    %br = linspace(-1, 5)';                  % Range to evaluate the beta posterior on
    %b1 = linspace(post.mu2(1) - sqrt(post.sigma2(1,1))*2.5, post.mu2(1) + sqrt(post.sigma2(1,1))*2.0);
    %b2 = linspace(post.mu2(2) - sqrt(post.sigma2(2,2))*2.5, post.mu2(2) + sqrt(post.sigma2(2,2))*2.0);
    b1 = linspace(1.3, 1.7, 500)';
    b2 = linspace(3.4, 3.6, 500)';
    [brx, bry] = meshgrid(b1,b2);
    br = [brx(:), bry(:)];
    bpost = post.pbeta(br);

    %plot(br, bpost);               % Plot beta posterior
    l = linspace(5, 15, 1e3);
    lpost = post.plambda(l);

    f = linspace(0, 0.95,1e3);
    fpost = post.pphi(f);
    end


 
    % 100 samples from a distribution different than the posterior
    bad_sample = [randn(100,2), exp(randn(100,1)), rand(100,1)*0.95];

    % 100 samples from same distribution as posterior

    % This would be a less manual way to draw from the posterior
    good_sample = draw_posterior_sample(param, post, 100); 


    % Draw from DRAM to verify
    param.nsimu = 100000;
    tic
    sol = dram_from_linver(param, false);
    toc

    good_sample = sol.chain(20001:500:end,:);


    % Draw sample from DRAM with bug added to sum-of-squares function

    solbad = dram_from_linver(param, true);
    
    bad_sample = solbad.chain(20001:500:end, :);

    % Confidence level
    alpha = 0.10;
    % Number of tests to run
    numtests = 50;

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

    save('dram_example4.mat', 'result_good', 'result_bad', 'good_sample', 'bad_sample', 'sol', 'solbad', 'G', 'y', 'e', 'post', '-mat')
