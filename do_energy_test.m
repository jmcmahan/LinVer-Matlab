function result = do_energy_test(versample, param, post, alpha, numtests, N1)
% result = do_energy_test(versample, param, post) : Do the energy statistic
% test 100 times for 100 samples from the posterior compared with the sample 
% provided in versample at confidence 0.05. Result is a structure:
%
% result.fail_ratio = ratio of tests which failed
% result.exact_sample = sample drawn from the exact posterior
%
% result = do_energy_test(versample, param, post, alpha) : Same as above but set
% confidence to alpha.
%
% result = do_energy_test(versample, param, post, alpha, numtests) : Same as above
% but do "numtests" number of tests. 
%
% result = do_energy_test(versample, param, post, alpha, numtests, N1) : Same
% as above but use N1 samples from the exact posterior.

if nargin < 4
    alpha = 0.05;
end

if nargin < 5
    numtests = 100;
end

if nargin < 6
    N1 = 100;
    N2 = length(versample);
end



s2 = versample;


if ~iscolumn(s2)
    s2 = s2';
end

numfail = 0;
disp('Starting energy tests...')
for j = 1:numtests
    % Drawn from the exact posterior
    s1 = draw_posterior_sample(param, post, N1);
    result.exactsamp((1 + (j-1)*N1):(j*N1), :) = s1;
    r = energy_dist_test(s1', s2', alpha);
    if r
        numfail = numfail + 1;
    end
    disp(sprintf('Test %d: %d', j, r));
end
result.fail_ratio = numfail / numtests;

end
