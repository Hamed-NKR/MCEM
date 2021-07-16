function [pp1, pp2] = CONNECT(pp1, pp2)
% "CONNECT" joins two approaching aggregates toward the minimum distance...
%     ...between their primaries.
% ----------------------------------------------------------------------- %
% 
% Inputs/Outputs:
%     pp1: Primary particle info array for aggregate 1
%     pp2: ~ for aggregate 2
% ----------------------------------------------------------------------- %

n_pp1 = size(pp1, 1); % Number of primary particles within aggregate 1
n_pp2 = size(pp2, 1); % ~ aggregate 2
pp1_temp = repelem(pp1, n_pp2, 1);
pp2_temp = repmat(pp2, n_pp1, 1);

dist_pp = sqrt(sum((pp1_temp(:,3:5) - pp2_temp(:,3:5)).^2, 2)) -...
    (pp1_temp(:,2) + pp2_temp(:,2)) ./ 2; % List of distances of primary...
        % ...particle centers
ind_cnt = find(dist_pp == min(dist_pp), 1); % Finding pairs with minimum...
    % ...distance
ind_cnt = [ind_cnt, pp1_temp(ind_cnt,1), pp2_temp(ind_cnt,1)]; % Adding...
    % ...the index numbers of the nearest pair in the overlapping matrix
ind_cnt = [ind_cnt, find(pp1(:,1) == ind_cnt(2)),...
    find(pp2(:,1) == ind_cnt(3))]; % Adding the actual indices of primaries

% Moving the aggregates to the point-touching status
dr = pp2(ind_cnt(5), 3:5) - pp1(ind_cnt(4), 3:5);
norm_dr = sqrt(sum(dr.^2, 2));
dr1 = dr ./ 2 - pp1(ind_cnt(4), 2) .* dr ./ (2 .* norm_dr);
dr2 = - dr ./ 2 + pp2(ind_cnt(5), 2) .* dr ./ (2 .* norm_dr);
pp1(:, 3:5) = pp1(:, 3:5) + repmat(dr1, n_pp1, 1);
pp2(:, 3:5) = pp2(:, 3:5) + repmat(dr2, n_pp2, 1);

end

