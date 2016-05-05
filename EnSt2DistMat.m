function e = EnSt2DistMat(ind1, ind2, dmat)
% Compute the energy distance statistic using the pre-computed distance
% matrix. The "ind1" and "ind2" should be permutations of indices of the
% pooled sample, each of the size of the two samples being tested. That is,
% X has some size N1, Y has some size N2, the pooled sample has some size
% N = N1 + N2, so ind1 is some permutation of 1:N of size N1, and ind2 are
% the remaining indices not chosen. 

n1 = length(ind1);
n2 = length(ind2);

N = n1 + n2;

diff1 = 0;
diff2 = 0;
diff3 = 0;

for i = 1:n1
    for j = 1:n2
        diff1 = diff1 + dist_from_dmat(ind1(i), ind2(j), dmat, N);
    end
end

for i = 1:n1
    for j = 1:n1
        diff2 = diff2 + dist_from_dmat(ind1(i), ind1(j), dmat, N);
    end
end

for i = 1:n2
    for j = 1:n2
        diff3 = diff3 + dist_from_dmat(ind2(i), ind2(j), dmat, N);
    end
end

e = (diff1*2/(n1*n2) - diff2/n1^2 - diff3/n2^2) * n1*n2 / N;