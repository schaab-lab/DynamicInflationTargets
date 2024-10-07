function p = define_parameters(varargin)


%% CALIBRATION

p.sigma = 0; % EIS set to 0 for analytical tractability (see paper)

% Parameters taken from Gali (2015):
p.beta = 0.99; % discount factor
p.phi = 5; % inverse Frisch elasticity
p.xi = 0.75; % Calvo parameter
p.epsilon = 9; % elasticity of substitution between goods

p.kappa = (1-p.xi) * (1-p.xi*p.beta) / p.xi * (p.sigma + p.phi); 

p.alpha = p.kappa / p.epsilon; % 0.75
p.lambda = (1-p.xi) * (1-p.xi*p.beta) / p.xi / (p.epsilon^2); % 0.05;

% Parameters taken from Ascari and Ropele (2007):
p.sigma_AR2007 = 1; % never used
p.epsilon_AR2007 = 11;

% Default shock persistence:
p.rho = 0.60;

% Additional parameters for r* application:
% p.lambda0 = 0.001;
% p.lambda1 = 0.01;
% p.rstar = 0.02;
% p.rho_bar = 0;
p.v0 = 0.1;
p.v1 = 0.4;


%% IRFS

% length of time grid:
p.T = 50;


%% COMPOSITE PARAMETERS
p.alpha_hat = p.alpha / (p.kappa^2);
p.lambda_hat = p.lambda / p.kappa;

% Section 4.1:
p.delta1 = 1/(2*p.alpha_hat*p.beta) * ( ...
           1 + p.alpha_hat*(1+p.beta) + p.v1 - ...
           sqrt((1+p.alpha_hat*(1+p.beta) + p.v1)^2 - 4*p.alpha_hat^2*p.beta));

p.delta2(1) = 1/p.alpha_hat * (1+p.alpha_hat*(1-p.beta*p.rho)) ...
               / (1-p.beta*p.delta1*p.rho) * p.delta1*p.v1;
p.delta2(2) = 1/p.alpha_hat * (1+p.alpha_hat*(1-p.beta*0)) ...
               / (1-p.beta*p.delta1*0) * p.delta1*p.v1;

p.delta0 = -p.delta1 * (1+p.alpha_hat*(1-p.beta)) ...
           / (p.alpha_hat * (1-p.beta*p.delta1)) * p.v0;
p.zeta = p.alpha_hat * p.beta / ((p.alpha_hat*p.beta + p.v1)*(1+p.alpha_hat) ...
         - p.alpha_hat^2*p.beta);

% Appendix B.2:
p.gamma2 = 1/(2*p.alpha_hat*p.beta) * ( ...
         1 + p.alpha_hat*(1+p.beta) - ...
         sqrt((1+p.alpha_hat*(1+p.beta))^2 - 4*p.alpha_hat^2*p.beta) );

p.gamma1(1) = p.gamma2 / (1-p.gamma2*p.beta*p.rho) ...
            * (p.rho-(1+p.alpha_hat)/(p.alpha_hat*p.beta));
p.gamma1(2) = p.gamma2 / (1-p.gamma2*p.beta*0) ...
            * (0-(1+p.alpha_hat)/(p.alpha_hat*p.beta));

%gamma1(1) = - p.alpha_hat/p.beta * p.gamma2 / (1-p.beta*p.rho*p.gamma2);
%gamma1(2) = - p.alpha_hat/p.beta * p.gamma2 / (1-p.beta*0*p.gamma2);
p.gamma0 = p.lambda_hat/p.alpha_hat * p.gamma2 / (1-p.beta*p.gamma2);



end

