function d = dist_from_dmat(i, j, dmat, N)
% Grab the distance from distance matrix "dmat", which is a vector storing
% the triangular elements of the pairwise distance matrix.

% Always have i < j
if j > i
    k = i;
    i = j;
    j = k;
end


% This formula just converts the i,j index into the vector format we've
% used. 

d = dmat(i + (j-1)*N - (j-1)*j/2);
end
