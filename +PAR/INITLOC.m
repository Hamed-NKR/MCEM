function pp_r = INITLOC(dom_size,pp_d)
% "INITLOC" randomly distributes the primary particles throughout the
% computational domain.

% Inputs are domain size and primary particle size array (mean+std).

n_pp = size(pp_d,1); % Number of primaries

% Initialization of the location array
pp_r = rand(n_pp,3) .* (repmat((dom_size)',n_pp,1) - repmat(pp_d,1,3)) +...
    (repmat(pp_d,1,3) ./ 2);

% Generating the "OVR" inputs:
ind_pps = (1:n_pp)'; % particle indices
ind_pps = [repelem(ind_pps,n_pp,1), repmat(ind_pps,n_pp,1)];
d_pps = [repelem(pp_d,n_pp,1), repmat(pp_d,n_pp,1)]; % size input
r_pps = [repelem(pp_r,n_pp,1), repmat(pp_r,n_pp,1)]; % location input
rmv = find(ind_pps(:,1) == ind_pps(:,2)); % removing self-identical pairs
ind_pps(rmv,:) = [];
d_pps(rmv,:) = [];
r_pps(rmv,:) = [];

ovrs = COL.OVR(r_pps, d_pps); % Checking initial overlapping...
% ...between the primaries
ind_err = 0; % Initializing error generation index

% Reinitializing overlapped particles
while ~ isempty(find(ovrs == 1, 1))
    ind_err = ind_err + 1; % Updating error index
    
    % Updating the location of overlapped pars
    ind_updt = ind_pps(ovrs == 1, 1); % Indices of updated particles
    ind_updt = unique(ind_updt); % removing repeating indices
    pp_r(ind_updt) = rand(size(ind_updt,1),3) .* (repmat((dom_size)',...
        size(ind_updt,1),1) - repmat(pp_d(ind_updt),1,3)) + ...
        (repmat(pp_d(ind_updt),1,3) ./ 2);
    r_pps = [repelem(pp_r,n_pp), repmat(pp_r,n_pp,1)];
    r_pps(rmv,:) = [];
    ovrs = COL.OVR(r_pps, d_pps); % Recheacking overlapping
    
    if ind_err > 10^3
        error("Error assigning random initial locations!\n")
    end
end

% n_pp = size(pp_d,1); % Number of primaries
% pp_r = zeros(n_pp,3); % Initialization of the location array
% ovr = 0; % overlap criterion
% i = 1; % primary particle loop index
% j = 0; % error generation index
% 
% % Finding random locations for the primaries
% while i <= n_pp
%     loc_tmp = (rand(3,1) .* (dom_size - pp_d(i))) + (pp_d(i) / 2); % a temporary location array
%     % checking the overlap with previously generated particles
%     if i > 1
%         k = 1; % Overlap checking loop index
%         while k <= i-1
%             ovr = COL.OVR([pp_r(k,:);loc_tmp'],[pp_d(k);pp_d(i)]);
%             if ovr == 1
%                 % overlap hapenning, the location needs recalculation.
%                 j = j+1;
%                 if j > 10^6
%                     error("Error assigning random initial locations!\n")
%                 else
%                     break
%                 end
%             else
%                 k = k+1;
%             end
%         end
%     end
%     if ovr == 1
%         continue % repeating the loop due to overlapping
%     else
%         pp_r(i,:) = loc_tmp'; % Assigning the unoverlapped location
%         j = 0 ; % Refreshing the error generation index
%         i = i+1; % moving forward to the next iteration
%     end
% end

end
