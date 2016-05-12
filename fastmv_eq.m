function y = fastmv_eq(phi, x)
% y = fastmv_eq(phi, x) - Return the matrix-vector product of x with the 
% Cholesky factorization of equicorrelation matrix 
% [ 1  phi phi   ...   phi]
% [phi  1  phi   ...   phi]
% [ .  .   .     ...    . ]
% [phi phi phi   ...    1 ]
%
% Note that if x is a sample from a standard normal distribution, y 
% is then equicorrelated Gaussian.
% 
% The product is computed using a recursive algorithm in O(N) time as opposed
% to the O(N^2) of the standard Cholesky algorithm. This can be used to make
% larger samples more feasible. 





y = zeros(size(x));

d = 1; o = phi;

% Running sum
s = 0;


for j = 1:length(x)
    y(j) = d*x(j) + s;
    s = s + o*x(j);
    dn = sqrt(d^2 - o^2); 
    o = (d - o)*o / dn;
    d = dn;
end


end
