function e = eval_noise(param)
% e = eval_noise(param) - Generate a set of random samples drawn from the
% Gaussian distribution specified by param


corrfunc = param.corrfunc;
N = param.N;
x = randn(N,1);
if N < 5e2;
    R = eval_corrfunc(param);
    L = chol(R)';
    e = L*x; 
else
    % Placeholder for a more efficient version 
    disp('N is larger than allowed by default, not attempting.')
    return
end

end
