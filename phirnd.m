function [y, x, err, phix, phiu, cdfphi, idc, newvals, e] = phirnd(post, param, M, options)
% y = phirnd(post, param, M) - Draw M samples from the phi posterior defined by
% the posterior structure "post" and parameters "param". The "options" structure
% can have the fields:
% options.cdfsize = number of points to interpolate the phi posterior, used in
%                   determining a roughly uniformly samples posterior.
% options.descdfsize = number of points in the uniformly sampled posterior. This
%                      should be smaller than "cdfsize".
%
% [y, x] = phirnd(post, param, M) - Same as above, but the uniform random sample
% used to build y is passed in x.


% This is a very simplistic way to find phi - cdf pairs that are sampled roughly
% uniformly in the cdf. There are better ways to do this, but this is quick and
% works well enough for our purposes. 
if nargin < 4
    options.cdfsize = 1e4;
    options.descdfsize = 1e3;
else
    if options.descdfsize > options.cdfsize
        disp('Warning - options.descdfsize should be smaller than options.cdfsize.');
        disp('Swapping options.cdfsize and options.descdfsize.');
        cdfsize = options.descdfsize;
        options.descdfsize = options.cdfsize;
        options.cdfsize = cdfsize;
    end
end

cdfsize = options.cdfsize;
descdfsize = options.descdfsize;

phirange = param.phirange;


phi = linspace(phirange(1), phirange(end), cdfsize)';
pphi = post.pphi(phi);
cdfphi = cumtrapz(phi, pphi);


% Uniformly spaced CDF output
descdf = linspace(0,1,descdfsize)'; 

% The values of phi roughly evenly spaced in the output
phix = zeros(size(descdf));
% The cdf at these values
phiu = zeros(size(descdf));
err = zeros(size(descdf));
previ = 1;
for j = 2:(length(descdf) - 1)
    % Adjustment to only search parts of the array that haven't been ruled out
    % (takes advantage of monotonicity of the cdf). 
    idx = previ:cdfsize; 
    [e, i] = min(abs(cdfphi(idx) - descdf(j)));
    i = i + previ - 1; 
    err(j) = e; 
    phix(j) = phi(i); 
    phiu(j) = cdfphi(i);
    previ = i; 
end

phix(1) = phirange(1);
phix(end) = phirange(end);
phiu(1) = 0;
phiu(end) = 1;


% Take a few Newton-iterations to make our input-output pairs for the CDF closer
% to the uniform output sampling we're going for. 

newtonits = 3;

% Desired output values
des = descdf;
% Initial input values
p = phix;
% Difference between current and desired
f = interp1(phi, cdfphi, p, 'spline') - des;
% Density = derivative of CDF
df = post.pphi(p);

for j = 1:newtonits
    p = p - f./df;
    f = interp1(phi, cdfphi, p, 'spline') - des;
    df = post.pphi(p);
end

% This is the error from the desired CDF for the final result
e = f;

% These are the final output values
newvals = interp1(phi, cdfphi, p, 'spline');

% Want indices with improved result (i.e., smaller error) and to eliminate any
% nan's arising from division-by-zero
idc = find(abs(e) < abs(err) & ~isnan(newvals));

phix(idc) = p(idc);
phiu(idc) = newvals(idc);
err(idc) = e(idc);

x = rand(M,1);

%y = interp1(phiu, phix, x, 'spline');
y = interp1(phiu, phix, x, 'linear');

end
