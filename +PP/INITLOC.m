function pp_r = INITLOC(dom_size,pp_d)
% This function randomly distributes the primary particles throughout...
% the computational domain.

% Inputs are domain size and primary particle size array (mean+std).

n_pp = size(pp_d,1); % Number of primaries
pp_r = zeros(n_pp,3); % Initialization of the location array
ovr = 0; % overlap criterion
i = 1; % primary particle loop index
j = 0; % error generation index

% Finding random locations for the primaries
while i <= n_pp
    loc_tmp = (rand(3,1) .* (dom_size - pp_d(i))) + (pp_d(i) / 2); % a temporary location array
    % checking the overlap with previously generated particles
    if i > 1
        k = 1; % Overlap checking loop index
        while k <= i-1
            ovr = COL.OVR([pp_r(k,:);loc_tmp'],[pp_d(k);pp_d(i)]);
            if ovr == 1
                % overlap hapenning, the location needs recalculation.
                j = j+1;
                if j > 10^6
                    error("Error assigning random initial locations!\n")
                else
                    break
                end
            else
                k = k+1;
            end
        end
    end
    if ovr == 1
        continue % repeating the loop due to overlapping
    else
        pp_r(i,:) = loc_tmp'; % Assigning the unoverlapped location
        j = 0 ; % Refreshing the error generation index
        i = i+1; % moving forward to the next iteration
    end
end

end
