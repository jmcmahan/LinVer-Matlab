# LinVer
LinVer is a reference Matlab implementation of a verification framework for Bayesian inference algorithms outlined in the technical note ["User Guidelines and Best Practices for CASL VUQ Analysis Using Dakota"](http://www.osti.gov/scitech/biblio/1177043). It is
based on a linear regression problem for which analytical or semi-analytical solutions are known.It provides a 
rigorous means of testing output chains of Markov chain Monte Carlo (MCMC) algorithms used for 
Bayesian inference are distributed correctly
via an implementation of a hypothesis test for equal distributions based on the energy distance
statistic (see [this paper](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.226.377) for details). While the main goal of the code is as a reference for those interested in 
implementing the framework in other verification software, it is also usable, as-is, as a basic verification tool. Mathematical details of the framework can be found in Appendix A of the previously mentioned ["technical note"](http://www.osti.gov/scitech/biblio/1177043) and will be elaborated on in a forthcoming paper. 

Be aware that LinVer is in development and may contain some bugs. The implementation
of the energy statistic test is not completely verified. The calculation of the true posteriors is believed
to be correct, however, for the basic cases. 

## Basic Example

The files demo\_case1.m, demo\_case2.m, and demo\_case3.m provide examples of how to calculate the true posteriors for the three different cases of unknowns (Case 1: regression parameters unknown, Case 2: regression parameters and scale parameter unknown, Case 3: regression parameters, scale parameter, and correlation parameter unknown). The following bit of code illustrates the setup of a simple verification problem. Additional information can be found in the documentation (doc/linver\_manual.pdf). 

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

