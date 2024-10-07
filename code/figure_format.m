
% FIGURE FORMATTING
fontname = 'Times';

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultLineLineWidth', 3.5);

set(groot, 'defaultAxesXGrid', 'on')
set(groot, 'defaultAxesYGrid', 'on')
set(groot, 'defaultAxesFontSize', 12)
set(groot, 'defaultAxesFontName', fontname);
set(groot, 'defaultTextFontName', fontname);

color_palette = [65, 105, 250; ...
                 201, 66, 79; ...
                 039,169,65] ./ 255;

set(groot, 'defaultAxesColorOrder', color_palette);
