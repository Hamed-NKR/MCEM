function dg = GYRATION(pp, n_pp)
% "GYRATION" computes the gyration diameter of the aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     pp: Primary particle information cell array
%     n_pp: Number distribution of primaries
% ----------------------------------------------------------------------- %
% 
% Output:
%     dg: Gyration diameter
% ----------------------------------------------------------------------- %

n_agg = size(n_pp,1); % Total number of aggregates
dg = zeros(n_agg,1); % Initializing the size array
r_com = PAR.COM(pp, n_pp); % Center of mass coordinates of aggregates

for i = 1 : n_agg
    
    a2 = sum((repmat(r_com(i,:), n_pp(i), 1) - pp{i}(:,3:5)).^2, 2);
        % Distance between center of primaries and aggregate center of mass
    radg2 = 0.6 .* (pp{i}(:,2) ./ 2).^2; % Radius of gyration of...
        % ...primaries squared
    dg(i) = 2 .* sqrt(sum(pp{i}(:,2).^3 .* (a2 + radg2)) ./...
        sum(pp{i}(:,2).^3)); % Gyration diameter
    
end

end
