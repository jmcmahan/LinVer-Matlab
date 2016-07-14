function L = fastchol_eq(phi, N)
% L = fastchol_eq(phi, N) - Return the lower-triangular Cholesky factorization
% of the equicorrelation matrix 
% [ 1  phi phi   ...   phi]
% [phi  1  phi   ...   phi]
% [ .  .   .     ...    . ]
% [phi phi phi   ...    1 ]
%
% The matrix is computed using a recursive algorithm in O(N) time as opposed
% to the O(N^2) of the standard Cholesky algorithm. This can be used to make
% larger samples more feasible. 



L = zeros(N);

% The entries of the Cholesky factor for the equicorrelated matrix
% can be computed with a recursive formula

d = zeros(N, 1);        % diagonal entries
o = zeros(N, 1);        % off-diagonal entries (N not used)

d(1) = 1; o(1) = phi;

L(1,1) = d(1);
L(2:end, 1) = o(1);

for j = 2:N
    d(j) = sqrt(d(j-1)^2 - o(j-1)^2);
    o(j) = (d(j-1) - o(j-1)) * o(j-1) / d(j);
    L(j,j) = d(j);
    L(j+1:end, j) = o(j);
end


end
