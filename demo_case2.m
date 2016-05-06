N = 300;                % Number of data points
Nbeta = 1;              % Number of regression parameters

% True parameters. These are used to generate data and that
% data is used to infer these parameters. 
beta = randn(Nbeta, 1)*10;
lambda = 100 + rand*200;
phi = 0.5;

% Non-informative prior only needs the type
prior_noninformative.type = 'noninformative';   

% Gaussian prior needs the mean and covariance specified
prior_gaussian.type = 'gaussian';       % Gaussian prior
prior_gaussian.sigma0 = eye(Nbeta);     % Prior unscaled covariance
prior_gaussian.mu0 = beta*1.1;          % Prior mean

% This is the design matrix. 
G = randn(N, Nbeta);          
G(:,1) = 1;


param1.N = N;
param1.Nbeta = Nbeta;
param1.G = G;
param1.prior = prior_noninformative;
param1.beta = beta;
param1.lambda = lambda;
param1.phi = phi;

% Correlation type
param1.corrfunc = 'none';


% Which parameters are unknown
param1.unknowns = 'beta_lambda';

% Generate the Gaussian observation error
e1 = eval_noise(param1);


% This is the data to fit the parameters to. 
param1.y = G*beta + e1 / sqrt(lambda); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncorrelated error 
% Non-informative prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

post1 = eval_posterior(param1); 
disp('******************************************************')
disp('Beta,lambda unknown, uncorrelated noise, uniform prior')
disp('******************************************************')
post1 = eval_posterior(param1);
result1 = do_energy_test(param1, post1);
disp('Done.')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncorrelated error 
% Gaussian prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copy everything...
param2 = param1;
% ... but change to Gaussian prior
param2.prior = prior_gaussian;

post2 = eval_posterior(param2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equi-correlated error 
% Non-informative prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param3 = param1;
param3.corrfunc = 'equal';
e3 = eval_noise(param3);
param3.y = G*beta + e3 / sqrt(lambda);

post3 = eval_posterior(param3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equi-correlated error 
% Gaussian prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param4 = param3;
param4.prior = prior_gaussian;

post4 = eval_posterior(param4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AR(1) correlated error 
% Non-informative prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param5 = param1;
param5.corrfunc = 'ar';
e5 = eval_noise(param5);
param5.y = G*beta + e5 / sqrt(lambda);

post5 = eval_posterior(param5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AR(1) correlated error 
% Gaussian prior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param6 = param5;
param6.prior = prior_gaussian;

post6 = eval_posterior(param6);
