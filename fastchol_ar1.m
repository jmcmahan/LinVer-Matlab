function L = fastchol_ar1(phi, N)
% L = fastchol_ar1(phi, N) - Return the lower-triangular Cholesky factorization
% of the AR(1) matrix 
% [1         phi       phi^2      ...  phi^(N-1) ]
% [phi       1         phi        ...  phi^(N-2) ]
% [phi^2     phi       1          ...  phi^(N-3) ]
% [ .        .         .          ...     .      ]
% [ .        .         .          ...     .      ]
% [ .        .         .          ...     .      ]
% [phi^(N-1) phi^(N-2) phi^(N-3)  ...     1      ]
%
% The matrix is computed using a recursive algorithm in O(N) time as opposed
% to the O(N^2) of the standard Cholesky algorithm. This can be used to make
% larger samples more feasible. 



L = zeros(N);

% Column arising in a pattern in the Cholesky matrix
c = ( phi.^(0:(N-1)) )';

L(:, 1) = c; 

c = c * sqrt(1 - phi^2);
for j = 2:N
    L(j:end, j) = c(1:(end-(j-1)));
end


end
