figure(1)
h = gcf;
if ~all(h.Position == [0, 0, 600, 600])
    h.Position = [0, 0, 600, 600];
end
set(h, 'color', 'white')

plotlabel = 'on';
ploteb = 'on';

dpp_uc = linspace(5, 65, 1000);
da_uc = 100 * (dpp_uc / 17.8).^(1 / 0.35);

p11 = plot(da_uc, dpp_uc, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '-.',...
    'LineWidth', 2.5); % Plot universal correlation
hold on

p12 = scatter(1e9 * da1, 1e9 * dpp1, 25, [0.8500 0.3250 0.0980],...
    'filled', 'o'); % Plot monodisperse aggs

% if strcmp(plotlabel, 'on')
%     text(1e9 * da1, 1e9 * dpp1, labels01, 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'right', 'FontSize', 8)
%     text(1e9 * da1, 1e9 * dpp1, labels02, 'VerticalAlignment', 'top',...
%         'HorizontalAlignment', 'left', 'FontSize', 8)
% end

p13 = scatter(1e9 * da21, 1e9 * dpp21, 30, [0.9290 0.6940 0.1250],...
    'filled', '^'); % Plot hybrid aggs - stage 1

std11 = PAR.MEANPP(pars1.pp);
labels11 = num2str(std11(:,2));
std12 = PAR.GEOMEANPP(pars1.pp);
labels12 = num2str(std12(:,2), '%.2f');

if strcmp(plotlabel, 'on')
%     text(1e9 * da21, 1e9 * dpp21, labels11, 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'right', 'FontSize', 8)
    text(1e9 * da21, 1e9 * dpp21, labels12, 'VerticalAlignment', 'baseline',...
        'HorizontalAlignment', 'left', 'FontSize', 8)
end

if strcmp(ploteb, 'on')
    e_p1 = dpp21 .* abs(std12(:,2) - 1);
    e_n1 = dpp21 .* abs(1 - 1 ./ std12(:,2));
    eb1 = errorbar(1e9 * da21, 1e9 * dpp21, 1e9 * e_n1, 1e9 * e_p1, '.');
    eb1.Color = [0.9290 0.6940 0.1250];
end

p14 = scatter(1e9 * da22, 1e9 * dpp22, 35, [0.4940 0.1840 0.5560],...
    'filled', 's'); % Plot hybrid aggs - stage 2

std21 = PAR.MEANPP(pars2.pp);
labels21 = num2str(std21(:,2));
std22 = PAR.GEOMEANPP(pars2.pp);
labels22 = num2str(std22(:,2), '%.2f');

if strcmp(plotlabel, 'on')
%     text(1e9 * da22, 1e9 * dpp22, labels21, 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'right', 'FontSize', 8)
    text(1e9 * da22, 1e9 * dpp22, labels22, 'VerticalAlignment', 'baseline',...
        'HorizontalAlignment', 'left', 'FontSize', 8)
end

if strcmp(ploteb, 'on')
    e_p2 = dpp22 .* abs(std22(:,2) - 1);
    e_n2 = dpp22 .* abs(1 - 1 ./ std22(:,2));
    eb2 = errorbar(1e9 * da22, 1e9 * dpp22, 1e9 * e_n2, 1e9 * e_p2, '.');
    eb2.Color = [0.4940 0.1840 0.5560];
end

p15 = scatter(1e9 * da23, 1e9 * dpp23, 30, [0 0.4470 0.7410],...
    'filled', 'v'); % Plot hybrid aggs - stage 3

std31 = PAR.MEANPP(pars3.pp);
labels31 = num2str(std31(:,2));
std32 = PAR.GEOMEANPP(pars3.pp);
labels32 = num2str(std32(:,2), '%.2f');

if strcmp(plotlabel, 'on')
%     text(1e9 * da23, 1e9 * dpp23, labels31, 'VerticalAlignment', 'bottom',...
%         'HorizontalAlignment', 'right', 'FontSize', 8)
    text(1e9 * da23, 1e9 * dpp23, labels32, 'VerticalAlignment', 'baseline',...
        'HorizontalAlignment', 'left', 'FontSize', 8)
end

if strcmp(ploteb, 'on')
    e_p3 = dpp23 .* abs(std32(:,2) - 1);
    e_n3 = dpp23 .* abs(1 - 1 ./ std32(:,2));
    eb3 = errorbar(1e9 * da23, 1e9 * dpp23, 1e9 * e_n3, 1e9 * e_p3, '.');
    eb3.Color = [0 0.4470 0.7410];
end

% axis equal
box on
set(gca, 'FontName', 'SansSerif', 'FontSize', 12, 'TickLength', [0.02 0.02])
xlabel({'\fontsize{14}d_a (nm)', '\fontsize{4} '},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
ylabel({'\fontsize{4} ', '\fontsize{14}d_p (nm)'},'interpreter','tex',...
    'FontName', 'SansSerif', 'FontWeight', 'bold')
ylim([5, 65])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend([p11, p12, p13, p14, p15], {'Universal correlation',...
    'Monodisperse', 'Hybrid - Young', 'Hybrid - Mid-age', 'Hybrid - Old'},...
    'Location', 'northwest', 'FontName', 'SansSerif', 'FontSize', 12);
title('Primary particle size vs projected area equivalent size',...
    'FontName', 'SansSerif', 'FontWeight', 'bold', 'FontSize', 16)
