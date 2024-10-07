function ll = figure_irf_panel(x, y, title_id, title_string, plot_xlabel, xlim_max, param)

N = 2;
ll = [];

%% OPTIONS
color_palette = [65, 105, 250; ...
                 201, 66, 79; ...
                 039,169,65] ./ 255;


% 31, 0, 168
% 91, 204, 10
% 255, 187, 41
% 8, 156, 94
% 69, 173, 87
% 87, 201, 66

% Combination 1:
% 63,23,189
% 170,38,0

% Combination 2:
% 66, 88, 201
% 87, 201, 66

% Combination 3:
% 41, 72, 186  42, 39, 179   65, 105, 250   31, 54, 140
% 201, 66, 79

% 3-color combo #1:
% 31, 54, 140
% 201, 66, 79
% 039,169,65   41, 161, 14

% combo #2:
% 65, 105, 250
% 201, 66, 79
% 110, 189, 60

% combo #3:
% 31, 54, 140
% 201, 66, 79
% 235, 201, 9

y_ticks = 3;
x_ticks = 6;

tick_font_size = 12;
x_label_font_size = 15;

line_width = 4;
% xlim_max = 40;

title_id_size = 17;
title_string_size = 17;


%% PLOT PANEL
hold on;
ll(1) = plot(x{1}, y{1}, ...
    'LineWidth', line_width, 'Color', color_palette(1, :), ...
    'LineStyle', '-');
ll(2) = plot(x{2}, y{2}, ...
    'LineWidth', line_width-2, 'Color', color_palette(2, :), ...
    'LineStyle', '--');
hold off; xlim([0, xlim_max]);
ax = gca; ax.TitleFontWeight = 'normal';
ylim = [ax.YLim(1), ax.YLim(2)];

a = get(gca, 'XTickLabel');
set(gca, 'FontSize', tick_font_size);
a = get(gca, 'YTickLabel');
set(gca, 'FontSize', tick_font_size);

ax.YAxis.Exponent=0;
%if max(abs(ylim)) >= 0.1 && max(abs(ylim)) < 1, ytickformat('%.1f'); end
%if max(abs(ylim)) >= 0.01 && max(abs(ylim)) < 0.1, ytickformat('%.2f'); end

set(gca, 'Xlim', [0, xlim_max]);
set(gca, 'Ylim', [ylim(1), ylim(2)]);

y_tick_step_size = (ax.YLim(2) - ax.YLim(1)) / (y_ticks-1);
set(gca, 'YTick', [ax.YLim(1) : y_tick_step_size : ax.YLim(2)])

x_tick_step_size = (ax.XLim(2) - ax.XLim(1)) / (x_ticks-1);
set(gca, 'XTick', [ax.XLim(1) : x_tick_step_size : ax.XLim(2)]);

title(['({\bf\fontsize{', num2str(title_id_size), '}', title_id, ...
       '})  {\fontsize{', num2str(title_string_size), '}', title_string, '}'], 'FontSize', 18);
if plot_xlabel, xlabel('Quarters', 'FontSize', x_label_font_size); end



%% FIGURE 2
%{
color_palette = [000,076,153; ...
                 255,188,000; ...
                 039,169,65; ...
                 000,205,108; ...
                 064,173,090; ...
                 202,091,035] ./ 255;

y_ticks = 5;
x_ticks = 6;

tick_font_size = 15;
x_label_font_size = 15;

line_width = 2;
xlim_max = min(15, param.T);

[~, id] = min(abs(param.t - xlim_max));
scatter_ticks = 10;
scatter_idx = 1 : round(id / scatter_ticks) : id;
[~, id] = min(abs(sim_RA.t - xlim_max));
scatter_ticks = 10;
scatter_idx2 = 1 : round(id / scatter_ticks) : id;

marker_size = 200;
% scatter(x,y,sz,'MarkerEdgeColor',[0 .5 .5],...
%               'MarkerFaceColor',[0 .7 .7],...
%               'LineWidth',1.5)

figure;
hold on;
l1 = plot(param.t, 100 * sim{1}.piw, 'LineWidth', line_width);
l2 = plot(sim_RA.t, 100 * sim_RA.piw, 'LineWidth', line_width);
scatter(param.t(scatter_idx), 100 * sim{1}.piw(scatter_idx), marker_size, ...
    'MarkerEdgeColor', color_palette(1, :), ...
    'MarkerFaceColor', color_palette(1, :), ...
    'LineWidth', 0.5);
scatter(sim_RA.t(scatter_idx2), 100 * sim_RA.piw(scatter_idx2), marker_size, ...
    'MarkerEdgeColor', color_palette(2, :), ...
    'MarkerFaceColor', color_palette(2, :), ...
    'LineWidth', 3, 'Marker', 'x');
hold off; xlim([0, xlim_max]);
ax = gca; ax.TitleFontWeight = 'normal';

y_tick_step_size = (ax.YLim(2) - ax.YLim(1)) / (y_ticks-1);
set(gca, 'YTick', [ax.YLim(1) : y_tick_step_size : ax.YLim(2)])

x_tick_step_size = (ax.XLim(2) - ax.XLim(1)) / (x_ticks-1);
set(gca, 'XTick', [ax.XLim(1) : x_tick_step_size : ax.XLim(2)])

a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'FontSize', tick_font_size);
a = get(gca, 'YTickLabel');
set(gca, 'YTickLabel', a, 'FontSize', tick_font_size);

title('this text is {\bfbold} this text is {\ititalizized}')
title('({\bf\fontsize{25}a}) {\fontsize{18}Optimal Interest Rate}', 'FontSize', 22);
xlabel('Quarters', 'FontSize', x_label_font_size); 
%}

end
