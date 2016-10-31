function e = eval_noise(param, scaled)
% e = eval_noise(param) - Generate a set of random samples drawn from the
% Gaussian distribution specified by param. NOTE: This is unscaled noise,
% meaning that it doesn't use param.lambda. If you are using this to
% generate observation error for the framework, then divide e by 
% sqrt(param.lambda). 
%
% e = eval_noise(param, scaled) - Same as above if scaled = false, 
% scale e by 1 / sqrt(param.lambda) if scaled = true.

if nargin < 2
    scaled = false;
end
    
corrfunc = param.corrfunc;
N = param.N;
x = randn(N,1);
if N <= 5e2;
    R = eval_corrfunc(param);
    if scaled
        L = chol(R)' / sqrt(param.lambda);
    else
        L = chol(R)'; 
    e = L*x; 
    end
else
    % Placeholder for a more efficient version 
    disp('N is larger than allowed by default, not attempting.')
    disp('Can change me in "eval_noise.m".');
    return
end

end
