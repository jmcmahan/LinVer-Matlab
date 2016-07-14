function [dmat] = dist_matrix(X, Y)
% Compute the distance matrix for pooled samples X and Y for use
% with the energy distance code. dmat is a vector holding the 
% triangular portions of the distance matrix. 
% X is d by n1, Y is d by n2, with d being the dimension, n1,n2 being
% the number of samples.

[d, N1] = size(X);
N2 = size(Y,2);

N = N1 + N2;

dmat = zeros((N+1)*N/2, 1);

Z = [X, Y];
k = 1;
for i = 1:N
    for j=i:N
        dmat(k) = norm(Z(:,i) - Z(:,j));
        k = k + 1;
    end
end
end