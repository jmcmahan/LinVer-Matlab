function R = eval_corrfunc(param, phi)
% R = eval_corrfunc(param) - Return the correlation matrix specified by the
% structure param (i.e., use the true phi value given by param.phi)
% 
% R = eval_corrfunc(param, phi) - Use the argument phi to determine the 
% correlation matrix R.

if nargin < 2
    phi = param.phi;
end

corrfunc = param.corrfunc;
N = param.N;

if strcmp(corrfunc, 'none')
    R = 1;  % This'll do the right thing in subsequent equations
elseif strcmp(corrfunc, 'equal')
    R = eye(N)*(1 - phi) + ones(N)*phi;
elseif strcmp(corrfunc, 'ar')
    R = zeros(N);
    for j = 1:N
        s = 1-j;
        R(:,j) = phi.^abs((s:(s+N-1)))';
    end
else
    disp('Error in eval_corrfunc: Unsupported correlation function.')
end

end
