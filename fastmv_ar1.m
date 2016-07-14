function y = fastmv_ar1(phi, x)
% y = fastmv_ar1(phi, x) - Return the matrix-vector product of x with the
% the Cholesky factor of the AR(1) matrix 
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



y = zeros(size(x));

y(1) = x(1);

for j = 2:length(x)
    y(j) = phi*y(j-1) + sqrt(1-phi^2)*x(j);
end

end
