function [y, x, cumphi] = phirnd_alt(post, param, M)
% y = phirnd(post, param, M) - Draw M samples from the phi posterior defined by
% the posterior structure "post" and parameters "param". 
%
% [y, x] = phirnd(post, param, M) - Same as above, but the uniform random sample
% used to build y is passed in x.


corrfunc = param.corrfunc;
xphi = post.phiparam.xphi;
wphi = post.phiparam.wphi;
pphi = post.phiparam.pphi;

phirange = param.phirange;

cumphi = zeros(length(wphi), 1);

for j = 1:length(cumphi)
    cumphi(j) = wphi(1:j)'*pphi(1:j);
end



% Find the closest cumulative-distribution value...
descdf = rand(M, 1);
D = repmat(descdf, 1, length(cumphi));
C = repmat(cumphi', M, 1);
[~, ind] = min(abs(D - C), [], 2);

% ...and approximate with linear interpolation to get the sample.
x1 = xphi(ind);

diff = descdf - cumphi(ind);
ip1 = min(max(ind + sign(diff), 1), length(descdf));
x2 = xphi(ip1);

deltac = abs(cumphi(ind) - cumphi(ip1));

% Binary search for the correct value
numits = 8; % number of binary iterations
y = x1;
for j = 1:numits
    y = y + sign(diff).*deltac/2; 
    diff = descdf - interp1(xphi, cumphi, y, 'spline');
    deltac = deltac/2;
end


x = descdf; 
end
