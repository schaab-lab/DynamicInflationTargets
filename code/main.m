%------------------------------------------------------------------------%
% 
% This code computes the impulse responses to lower bound, r*, and NKPC
% slope shocks under a dynamic inflation target.
% 
% Code written by Christopher Clayton and Andreas Schaab.
% Current version: October 2024. First version: September 2022.
% 
% Please cite our paper:
%   Clayton, C. and A. Schaab. A Theory of Dynamic Inflation Targets.
%   Conditionally accepted at American Economic Review. 
% Thanks!
% 
% This code uses the following function for an aesthetic axis break:
% Chris McComb (2024). Aesethetic Axis Breaks (https://github.com/cmccomb/break_axis), GitHub.
% 
%------------------------------------------------------------------------%

clear
close all
clc

diary ./output/output.log
diary on

figure_format;

fprintf('Running algorithm:\n')
run_time = tic;


%% PARAMETERS
p = define_parameters();


%% TIME GRID
t = 1:p.T; 


%% SECTION 4.1: r*
% v0 = 1/(2*p.beta) * (p.lambda0 + p.lambda1*p.rho_bar) ...
%       / (p.rho_bar + p.rstar);
% v1 = 1/(2*p.beta) * p.lambda1 / (p.rho_bar + p.rstar);

% Initialize:
[tau, nu, pi, theta] = deal(zeros(p.T, 2)); 

% Shock:
theta(1, :) = -0.05;
for n = 1:p.T-1
    theta(n+1, 1) = p.rho * theta(n);
    theta(n+1, 2) = 0;
end

% Dynamic inflation target: steady state
nu_RSS = p.delta0 / (1-p.delta1);
tau_RSS = p.zeta * ( ...
            p.v0 + p.delta0 + (p.delta1 - (p.alpha_hat*p.beta + p.v1) ...
            / (p.alpha_hat*p.beta)) * nu_RSS );

% Dynamic inflation target: transition dynamics
nu(1, :) = p.delta0 + p.delta2.*theta(1, :) + p.delta1 * nu_RSS;
for n = 2:p.T
    nu(n, :) = p.delta0 + p.delta2.*theta(n, :) + p.delta1*nu(n-1, :);
end

b = - nu;

tau(:, 1) = p.zeta * ( ...
                p.v0 + p.delta0 ...
                + (p.delta1 - (p.alpha_hat*p.beta + p.v1) ...
                  / (p.alpha_hat*p.beta)) * nu(:, 1) ...
                + (p.delta2(:, 1) - p.v1) * p.rho * theta(:, 1) );
tau(:, 2) = p.zeta * ( ...
                p.v0 + p.delta0 ...
                + (p.delta1 - (p.alpha_hat*p.beta + p.v1) ...
                  / (p.alpha_hat*p.beta)) * nu(:, 2) );

pi(1, :) = p.zeta * (nu(1, :) ...
    - (p.alpha_hat*p.beta+p.v1)/(p.alpha_hat*p.beta) * nu_RSS ...
    + p.v0 - p.v1 * theta(1, :) );
for n = 2:p.T
    pi(n, :) = p.zeta * (nu(n, :) ...
        - (p.alpha_hat*p.beta+p.v1)/(p.alpha_hat*p.beta) * nu(n-1, :) ...
        + p.v0 - p.v1 * theta(n, :) );
end

% Plot figure:
xx{1} = t; xx{2} = t; x_max = 20;
figure;
subplot(2, 2, 1);
y{1} = b(:, 1); y{2} = b(:, 2);
title_string = 'Target Flexibility';
ll = figure_irf_panel(xx, y, 'A', title_string, 0, x_max, p);

subplot(2, 2, 2);
y{1} = tau(:, 1); y{2} = tau(:, 2);
title_string = 'Target Level';
figure_irf_panel(xx, y, 'B', title_string, 0, x_max, p);

subplot(2, 2, 3);
y{1} = pi(:, 1); y{2} = pi(:, 2);
title_string = 'Inflation';
figure_irf_panel(xx, y, 'C', title_string, 0, x_max, p);

subplot(2, 2, 4);
y{1} = theta(:, 1); y{2} = theta(:, 2);
title_string = 'Shock (r*)';
figure_irf_panel(xx, y, 'D', title_string, 0, x_max, p);

legend(ll, {'{\bfPersistent}', '{\bfTransitory}'}, ...
    'Location', 'NorthEast', 'box', 'off', 'FontSize', 14);

set(gcf,'Position',[500, 500, 1000, 350]);
print('./output/rstar.eps', '-depsc', '-painters', '-noui', '-r600');


%% SECTION 4.2: SLOPE NKPC

% Initialize:
[pi, theta] = deal(zeros(p.T, 2));

% Shocks:
theta(1, :) = 1.05;
for n = 1:p.T-1
    theta(n+1, 1) = 1 - p.rho + p.rho * theta(n);
    theta(n+1, 2) = 1;
end

% Dynamic inflation target:
nu = theta / p.kappa;
b = - nu;

tau = (1-[p.rho 0]) .* (1/p.kappa - theta/p.kappa);

for n = p.T:-1:2
    pi(n, :) = theta(n, :)/p.kappa - theta(n-1, :)/p.kappa;
end
pi(1, :) = theta(1, :)/p.kappa - 1/p.kappa;

% Plot figure:
xx{1} = t; xx{2} = t; x_max = 20;
figure;
subplot(2, 2, 1);
y{1} = b(:, 1); y{2} = b(:, 2);
title_string = 'Target Flexibility';
ll = figure_irf_panel(xx, y, 'A', title_string, 0, x_max, p);

subplot(2, 2, 2);
y{1} = tau(:, 1); y{2} = tau(:, 2);
title_string = 'Target Level';
figure_irf_panel(xx, y, 'B', title_string, 0, x_max, p);

subplot(2, 2, 3);
y{1} = pi(:, 1); y{2} = pi(:, 2);
title_string = 'Inflation';
figure_irf_panel(xx, y, 'C', title_string, 0, x_max, p);

subplot(2, 2, 4);
y{1} = p.kappa ./ theta(:, 1); y{2} = p.kappa ./ theta(:, 2);
title_string = 'Shock (NKPC Slope)';
figure_irf_panel(xx, y, 'D', title_string, 0, x_max, p);

legend(ll, {'{\bfPersistent}', '{\bfTransitory}'}, ...
    'Location', 'SouthEast', 'box', 'off', 'FontSize', 14);

set(gcf,'Position',[500, 500, 1000, 350]);
print('./output/slope.eps', '-depsc', '-painters', '-noui', '-r600');


%% SECTION 5.1: COMMITMENT CURVE
k = 1:50;

% Trend inflation = 1%
gamma = 1.01;

beta_tilde = (gamma-1) * p.beta * (1-p.xi*gamma^(p.epsilon_AR2007-1)) ...
             * (p.epsilon_AR2007-1);
delta_tilde = p.xi * p.beta * gamma^(p.epsilon_AR2007-1);

beta_star = beta_tilde / (beta_tilde + p.beta*gamma);
delta_star = delta_tilde / p.beta;

% y1 = beta_star * delta_star.^(k-1);
y1 = beta_star*delta_star/(1-delta_star) / ...
    (1 + beta_star*delta_star/(1-delta_star)) * delta_star.^(k-2);
y1(1) = 1; % quasi-hyperbolic jump


% Trend inflation = 2%
gamma = 1.02;

beta_tilde = (gamma-1) * p.beta * (1-p.xi*gamma^(p.epsilon_AR2007-1)) ...
             * (p.epsilon_AR2007-1);
delta_tilde = p.xi * p.beta * gamma^(p.epsilon_AR2007-1);

beta_star = beta_tilde / (beta_tilde + p.beta*gamma);
delta_star = delta_tilde / p.beta;

% y1 = beta_star * delta_star.^(k-1);
y2 = beta_star*delta_star/(1-delta_star) / ...
    (1 + beta_star*delta_star/(1-delta_star)) * delta_star.^(k-2);
y2(1) = 1; % quasi-hyperbolic jump


% Plot figure:
y_ticks = 3;
x_ticks = 6;
tick_font_size = 12;
x_label_font_size = 15;
line_width = 4;
xlim_max = 15;
title_id_size = 17;
title_string_size = 17;

figure; hold on;
l1 = plot(k(2:end), y1(2:end), ...
    'LineWidth', line_width, 'Color', color_palette(1, :), ...
    'LineStyle', '-'); % xlim([1, xlim_max]);
l2 = plot(k(2:end), y2(2:end), ...
    'LineWidth', line_width, 'Color', color_palette(2, :), ...
    'LineStyle', '-'); % xlim([1, xlim_max]);
scatter(1, 0.25, 150, 'LineWidth', 3, 'MarkerEdgeColor', color_palette(1, :));
scatter(1, 0.25, 150, 'LineWidth', 3, 'MarkerEdgeColor', color_palette(2, :));
hold off;

break_axis('position', 0.18)
yticks([0 0.05 0.1 0.15 0.25])
yticklabels({'0', '0.05', '0.1', '0.15', '1'})

xlabel('$k$', 'Interpreter', 'Latex');

set(gcf,'Position',[500, 500, 550, 300]);

legend([l2, l1], {'$\bar \pi = 2\%$', '$\bar \pi = 1\%$'}, ...
    'Location', 'NorthEast', 'box', 'off', 'FontSize', 18, 'Interpreter', 'Latex');

print('./output/commitment_curve.eps', '-depsc', '-painters', '-noui', '-r600');


%% APPENDIX B.2: LOWER BOUND

% Initialize:
[nu, tau, theta] = deal(zeros(p.T, 2));

% Shock:
theta(1, :) = 0.05;
for n = 1:p.T-1
    theta(n+1, 1) = p.rho * theta(n);
    theta(n+1, 2) = 0;
end

% Dynamic inflation target: steady state
nu_RSS = p.gamma0 / (1-p.gamma2);
tau_RSS = 0;

% Dynamic inflation target: transition dynamics
nu(1, :) = p.gamma0 + p.gamma1.*theta(1, :) + p.gamma2 * nu_RSS;
for n = 2:p.T
    nu(n, :) = p.gamma0 + p.gamma1.*theta(n, :) + p.gamma2*nu(n-1, :);
end

for n = 1:p.T
    tau(n, :) = p.gamma0 + (p.gamma2-1)*nu(n,:) ...
                + (p.gamma1 + 1/p.beta).*[p.rho 0].*theta(n,:);
end

nulag = [nu_RSS nu_RSS];
nulag = [nulag; nu(1:end-1,:)];

pi = nu - nulag + 1/p.beta*theta;
b = - nu;

% Plot figure:
xx{1} = t; xx{2} = t; x_max = 20;
figure;
subplot(2, 2, 1);
y{1} = b(:, 1); y{2} = b(:, 2);
title_string = 'Target Flexibility';
ll = figure_irf_panel(xx, y, 'A', title_string, 0, x_max, p);

subplot(2, 2, 2);
y{1} = tau(:, 1); y{2} = tau(:, 2);
title_string = 'Target Level';
figure_irf_panel(xx, y, 'B', title_string, 0, x_max, p);

subplot(2, 2, 3);
y{1} = pi(:, 1); y{2} = pi(:, 2);
title_string = 'Inflation';
figure_irf_panel(xx, y, 'C', title_string, 0, x_max, p);

subplot(2, 2, 4);
y{1} = theta(:, 1); y{2} = theta(:, 2);
title_string = 'Shock (Lower Bound)';
figure_irf_panel(xx, y, 'D', title_string, 0, x_max, p);

legend(ll, {'{\bfPersistent}', '{\bfTransitory}'}, ...
    'Location', 'NorthEast', 'box', 'off', 'FontSize', 14);

set(gcf,'Position',[500, 500, 1000, 350]);
print('./output/lower_bound.eps', '-depsc', '-painters', '-noui', '-r600');


%% FINISH
run_time = toc(run_time); 
fprintf('\n\nAlgorithm converged. Run-time of: %.2f seconds.\n', run_time);

diary off
