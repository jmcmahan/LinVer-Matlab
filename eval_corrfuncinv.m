function R = eval_corrfuncinv(param, phi)
% R = eval_corrfuncinv(param) - Return the inverse correlation matrix 
% specified by the structure param (i.e., use the true phi value given by 
% param.phi)
% 
% R = eval_corrfuncinv(param, phi) - Use the argument phi to determine the 
% inverse correlation matrix R.

if nargin < 2
    phi = param.phi;
end

corrfunc = param.corrfunc;
N = param.N;

if strcmp(corrfunc, 'none')
    R = 1;
elseif strcmp(corrfunc, 'equal')
    R = (eye(N) - ones(N)*phi/(1 + (N-1)*phi)) / (1-phi);
elseif strcmp(corrfunc, 'ar')
    % Set according to the near-pattern...
    R = diag(ones(N,1)*(1+phi^2)) - phi*diag(ones(N-1,1), 1) ...
                                  - phi*diag(ones(N-1,1),-1);
    % ...then make corrections
    R(1,1) = 1; R(N,N) = 1; 
    R = R / (1-phi^2);

else
    disp('Error in eval_corrfunc: Unsupported correlation function.')
end

end
