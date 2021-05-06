function par = GROW(par)
% "GROW" monitors collision of particles and updates the particle structure
% based on the post-aggregation data.

% Input is the particle structure.

n_par = size(par.d,1);  % Total number of particles

for i = 1 : n_par - 1
    for j = i+1 : n_par
        % Checking collision occurrence
        if COL.OVR([par.r(i,:); par.r(j,:)], [par.d(i); par.d(j)])
                        
        end
    end
end

end

