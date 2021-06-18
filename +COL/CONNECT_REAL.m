function [pp1, pp2, colstat] = CONNECT_REAL(pp1, v_par1, pp2, v_par2)
% "CONNECT_REAL" sticks two approaching aggregates to each other based...
    % ...on their real physical status including their relative...
    % ...locations and velocities.
% ----------------------------------------------------------------------- %

% Inputs/Outputs:
    % pp"i": Primary particle info cell array for aggregate "i"
    % v_par"i": Center of mass velocity of aggregate "i"
    % stat: Returns the collision status of aggregate pairs 
% ----------------------------------------------------------------------- %

n_pp1 = size(pp1, 1); % Number of primary particles within aggregate 1
n_pp2 = size(pp2, 1); % ~ aggregate 2

% Preparing for vectorized operations
pp1_temp = repelem(pp1, n_pp2, 1);
v_par1_temp = repmat(v_par1, n_pp1 * n_pp2, 1);
pp2_temp = repmat(pp2, n_pp1, 1);
v_par2_temp = repmat(v_par2, n_pp1 * n_pp2, 1);

% Calculating the minimal distances and assocaited timesteps
del_r = COL.DYNDIST([pp1_temp(:,3:5), pp2_temp(:,3:5)],...
    [v_par1_temp, v_par2_temp]);
ind_cnt = find(del_r == min(del_r), 1); % Finding pairs with the...
    % ...lowest minimal distance

% Checking if there is a chance for the aggregates to collide sometime...
    % ...in the futre
if min(del_r) > (pp1_temp(ind_cnt,2) + pp2_temp(ind_cnt,2)) / 2
    colstat = 0; % Aggregates do not have a chance to collide
else
    colstat = 1; % Aggregates migh be able to collide
    ind_cnt = [ind_cnt, pp1_temp(ind_cnt,1), pp2_temp(ind_cnt,1)];
        % Adding the index numbers of the nearest pair
    ind_cnt = [ind_cnt, find(pp1(:,1) == ind_cnt(2)),...
        find(pp2(:,1) == ind_cnt(3))];
    
    dr_cnt = COL.PNTOUCH([pp1(ind_cnt(4), 3:5),...
        pp2(ind_cnt(5), 3:5)], [v_par1, v_par2], [pp1(ind_cnt(4), 2),...
        pp2(ind_cnt(5), 2)]); % The displacement needed for the nearest...
            % ...primaries to point-touch
    
    % Updating primary particle locations
    pp1(:, 3:5) = pp1(:, 3:5) + repmat(dr_cnt(1, 1:3), n_pp1, 1);
    pp2(:, 3:5) = pp2(:, 3:5) + repmat(dr_cnt(1, 4:6), n_pp2, 1);
end

end

