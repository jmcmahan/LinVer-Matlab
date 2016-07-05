function result = demo_energy_test(param, post, alpha, numtests, N1, N2)
% result = do_energy_test(param, post) : Do the energy statistic
% test 100 times for 100 samples from the posterior and 100 samples
% computed via MCMC at confidence 0.05. Result is a structure like so:
%
% result.exactchain = samples from the exact posterior
% result.mcmcsamp = sample from the MCMC posterior
% result.chain = entire MCMC chain
% result.accept_ratio = accept ratio of MCMC chain
% result.fail_ratio = ratio of tests which failed
%
% result = do_energy_test(param, post, alpha) : Same as above but set
% confidence to alpha.
%
% result = do_energy_test(param, post, alpha, numtests) : Same as above
% but do "numtests" number of tests. 
%
% result = do_energy_test(param, post, alpha, numtests, N1, N2) : Same
% as above but use N1 samples from the exact posterior and N2 samples
% from the MCMC chain.

if nargin < 3
    alpha = 0.05;
end

if nargin < 4
    numtests = 100;
end

if nargin < 6
    N1 = 100;
    N2 = 100;
end



% We'll take every 100th iterate as the MCMC sample
[chain, aratio] = do_simple_mcmc(param, post, 1 + N2*100); 

s2 = chain(2:100:end, :); 

result.chain = chain;
result.accept_ratio = aratio;
result.exactsamp = zeros(numtests*N1, size(s2,2));
result.mcmcsamp = s2;


res_etest = do_energy_test(s2, param, post, alpha, numtests, N1, N2)
result.exactsamp = res_etest.exactsamp;
result.fail_ratio = res_etest.fail_ratio;

end
