function [reject, pval, estat, B] = energy_dist_test(s1, s2, alpha, B)
% Do the two-sample energy statistic test for equality of distributions.
% s1 = d-by-N1 set of N1 samples of d-dimensional variables
% s2 = d-by-N2 set of N2 samples of d-dimensional variables
% alpha = significance level to use. The hypothesis that the distributions
%         are equal is rejected if more than 100*(1-alpha) percent of
%         the permutation samples from the pooled data are smaller than
%         the observed value.
% B  = number of permutation trials to use for approximating the 
%      distribution of the test statistic. If not provided, is
%      automatically computed from alpha. 
%
% The test compares the data in s1 to the data in s2 to see if the samples
% appear to come from the same distribution.
%
% reject = True if rejected at alpha significance level, False otherwise.
%          i.e., if true, the data supports the notion that s1 and s2 are
%          from different distributions at significance level alpha.
%          If false, it is undetermined.
% pval = ratio of permutation samples with energy distance greater than or 
%        equal to the observed energy distance (includes the observed).
% estat = vector of energy statistics for the B+1 permutation samples
%         (includes the non-permuted s1, s2)


if nargin < 4
    B = 5*ceil(1/alpha)-1;
end

if 1/(B+1) >= alpha
    disp('Warning: Choice of B may not be appropriate for desired signifcance, alpha');
end

% Note, need the same d for s1 and s2. 
[d, n1] = size(s1);
[d, n2] = size(s2);

N = n1+n2;

pool = [s1, s2];

estat = zeros(B+1, 1);
eobs = EnSt2(s1, s2);
estat(1) = eobs;

for j = 2:B+1
    idx = randsample(N, N);
    x = pool(:, idx(1:n1));
    y = pool(:, idx(n1+1:end));
    estat(j) = EnSt2(x, y);
end

pval = sum(eobs < estat) / (B+1);

if alpha > pval
    reject = true;
else
    reject = false;
end

end

