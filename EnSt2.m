function [e] = EnSt2(X, Y)
%EnSt2 Energy statistic for samples x and y
%   Computes the energy statistic for samples X and Y which are 
%   d by n1 and d by n2 matrices, respectively containing n1, n2 
%   observations of the d-dimensional variables.

[d, N1] = size(X);
N2 = size(Y,2);
n = N1 + N2;

% Implementation will always assume x has more samples
if N1 < N2
    x = Y;
    n1 = N2;
    y = X;
    n2 = N1;
else
    x = X;
    n1 = N1;
    y = Y;
    n2 = N2;
end


diff1 = 0;
diff2 = 0;
diff3 = 0;
% All pairwise subtractions implemented by circularly shifting 
for j = 1 : n1    
    xs = circshift(x, [0 j]);
    diff1 = diff1 + sumdists(xs(:,1:n2), y);
    diff2 = diff2 + sumdists(xs, x);
end

for j = 1 : n2
    ys = circshift(y, [0 j]);
    diff3 = diff3 + sumdists(ys, y);
end

e = (diff1*2/(n1*n2) - diff2/n1^2 - diff3/n2^2) * n1*n2 / n;

end



function v = sumdists(x, y)
% Computes the sum of the distances between the given sets of 
% observations. The number of observations (i.e., number of columns in
% this code) must be equal.
    v = sum(sqrt(sum(abs(x - y).^2,1)));
end
