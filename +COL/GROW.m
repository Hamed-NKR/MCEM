function [pp, agg] = GROW(pp,agg)
% "GROW" generates new aggregates upon collision of smaller aggregates or
% primaries.

% Inputs are primary particle and aggregate structures.

n_pp = size(pp.d,1);  % Number of primaries
n_agg = size(agg.d,1);  % Number of aggregates
n_par = n_pp + n_agg;  % Total number of particles
par_d = [pp.d; agg.d];
par_r = [pp.r; agg.r];

for i = 1 : n_par - 1
    for j = i+1 : n_par
        % Checking collision occurrence
        if COL.OVR([par_r(i,:); par_r(j,:)], [par_d(i); par_d(j)])
            if (i <= n_pp) && (j <= n_pp)
                % particle-particle aggregation
                pp = PPA(pp, ind_pp1, ind_pp2);
            elseif (i <= n_pp) && (j > n_pp)
                % particle-cluster aggregation
                [pp, agg] = PCA(pp, agg, ind_pp, ind_agg);
            else
                % cluster-cluster aggregation
                agg = CCA(agg, ind_agg1, ind_agg2);
            end
        end
    end
end

end

