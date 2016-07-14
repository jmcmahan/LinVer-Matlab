    param.N = 300;                          % Number of data points
    param.beta = [1.5 ];                    % True regression parameters
    param.Nbeta = length(param.beta);       % Number of regression parameters
    param.lambda = 10.0;                    % True scale parameter 
    param.phi = 0.5;                        % True correlation parameter
    param.corrfunc = 'equal';               % Correlation function 
    param.prior.type = 'noninformative';    % Use non-informative prior
    param.unknowns = 'beta';                % Compute posterior for the regression parameters only
    param.G = randn(param.N, param.Nbeta);  % Design matrix...
    param.G(:,1) = 1;                       % ...with first regression parameter a bias term.

    y0 = param.G * param.beta;              % Error-free observation data

    e = eval_noise(param);                  % Generate unscaled observation error
    param.y = y0 + e / sqrt(param.lambda);  % Add error to create data for calibration
    
    post = eval_posterior(param);           % Compute posterior

    br = linspace(-1, 5)';                  % Range to evaluate the beta posterior on
    plot(br, post.pbeta(br));               % Plot beta posterior

 
    % 100 samples from a distribution different than the posterior
    bad_sample = randn(100,1);

    % 100 samples from same distribution as posterior
    good_sample = post.mu2 + sqrt(post.sigma2 / param.lambda) * randn(100,1); 

    % This would be a less manual way to draw from the posterior
    %good_sample = draw_posterior_sample(param, post, 100); 

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
