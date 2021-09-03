function h_kin = PLOTKINETICS(parsdata)
% "PLOTKINETICS" displays the population-based kinetics of particle...
%	...collision and growth over time.
% ----------------------------------------------------------------------- %
%
% Input:
%   tbl: Table of particle properties
% ----------------------------------------------------------------------- %
%
% Output:
%   h_kin: The output figure handle
% ----------------------------------------------------------------------- %

% Initializing the plot outlines
hold off
clf('reset')
h_kin = gcf;
figure (h_kin);
if ~all(h_kin.Position == [0, 0, 2000, 892.1])
    h_kin.Position = [0, 0, 2000, 892.1]; % Setting the figure size
end
if ~all(h_kin.Color == [1 1 1])
    h_kin.Color = [1 1 1]; % Background color
end

% Changes of number and collision frequency over time
subplot(1,2,1)
yyaxis left
plot(parsdata.t, parsdata.ntot);
xlabel('t (s)')
ylabel('ntot (#)')
yyaxis right
plot(parsdata.t, parsdata.beta);
ylabel('beta (1/s)')
title ('Total number vs. time')
% set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

% Size distributions of particles mass over time
subplot(1,2,2)
n_data = length(parsdata.t); % Number of datasets to be plotted
% Setting a color distribution for the second plot lines
cl = jet;
ii1 = unique(round(10 + (length(cl) - 10) .* (0 : 1 / 10 : 1)'));
ii2 = unique(round(1 : (n_data - 1) / 10 : n_data)');
txt = cell(length(ii2),1);
for i = 1 : length(ii2)
    plot(parsdata.dm_dlogdv{ii2(i)}(:,2),...
        parsdata.dm_dlogdv{ii2(i)}(:,1),'Color', cl(ii1(i),:));
    txt{i} = "t = " + num2str(parsdata.t(ii2(i)), '%1.1e') + " (s)";
        % Time legend
    hold on
end
legend(txt)
set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
xlabel('dv (m)')
ylabel('dm/dlogdv (-)')
title ('Mass based size distribution')
hold off

if nargout == 0
    clear h_kin;  % Deleting figure handle if not requested as an output
end

end

