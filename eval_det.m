function d = eval_det(param, phi)
% d = eval_det(param, phi) - Evaluate the determinant for the correlation function 
% specified by % param at the correlation parameters given by the vector phi
% phi - vector of values at which to evaluate the determinant
% param - parameter structure containing a field "corrfunc". This can be one of
% three strings : "none", "equal", "ar", which determines the correlation
% function. It should also have the field "N" which specifies the size of the
% data. 



N = param.N;
corrfunc = param.corrfunc;
if strcmp(corrfunc, 'none')
    d = 1;
elseif strcmp(corrfunc, 'equal')
    d = ((1 - phi).^(N-1)) .* (1 + (N-1)*phi);
elseif strcmp(corrfunc, 'ar')
    d = (1 - phi.^2).^(N-1);
else
    disp('Error in eval_det: Unsupported correlation function.')
end

end
