function h = PCF_MULTI(g, r, lbl)
% PCF_MULTI gets multiple pair-correlation function (PCF) distributions...
%   ...from different aggregates/scenarios and visualizes them for...
%   ...comparison in a single PCF vs. radial location plot. 
% ----------------------------------------------------------------------- %
% 
% Inputs:
%   g: A cell array containing data of radial vector of averaged (over...
%       ...primary particles) PCF within and nearby different aggregates.
%   r: The raidal logarithmic points of calculation for PCF starting...
%       ...from inside a central primary particle. This can also be...
%       ...different for each aggregates.
%   lbl: A cell array of labels to be used for each aggregate data.
% ----------------------------------------------------------------------- %
% 
% Outputs:
%   h: PCF plot figure handle
% ----------------------------------------------------------------------- %

% initialize figure 
figure;
h = gcf;
h.Position = [0, 0, 700, 700];
set(h, 'color', 'white');

for i = 1 : length(g)
    if i < 5
        plot(r{i}, g{i});
    else
        plot(r{i}, g{i}, ':');
    end
    hold on
end

box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$\overline{r}$ (-)', 'interpreter', 'latex', 'FontSize', 20)
ylabel('$\overline{g}$($\overline{r}$) (-)', 'interpreter', 'latex', 'FontSize', 20)
legend(lbl, 'interpreter', 'latex', 'FontSize', 12, 'Location', 'northoutside',...
    'Orientation', 'horizontal', 'NumColumns', 2);

end

