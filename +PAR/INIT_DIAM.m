function [pp_d, pp_n] = INIT_DIAM(n_par, n_pp, d_pp)
% "INIT_DIAM" assigns a normal distribution to the global primary...
%     ...particle sizes and another normal size distribution to the...
%     ...primaries within the aggregates.
% ----------------------------------------------------------------------- %
% 
% Inputs:
%     n_par: The total number of (independent) particles.
%     n_pp: The number and standard deviation of primaries within each...
%         ...particle (a 2*1 vector).
%     d_pp: The global mean size and standard deviation of primaries...
%         ...as well as the size standard deviation within the aggregates.
% ----------------------------------------------------------------------- %
% 
% Outputs:
%     pp_d: Primary particle sizes
%     pp_n: Primary particle number distribution (numbers within each...
%         ...particle) 
% ----------------------------------------------------------------------- %

pp_d = zeros(n_par,1); % declaring the diameter array
pp_n = ones(n_par,1); % declaring the number array
chk = find(pp_d <= 0,1); % looping criterion
i = 0; % error generation index

% Generating a normally distributed series of random diameters and...
% ...internal primary numbers.
while ~ isempty(chk)
    pp_d = lognrnd(log(d_pp(1)), log(d_pp(2)), [n_par,1]); % Global size distribution (lognormal)
    chk = find(pp_d <= 0,1); % check all the diameters to be positive
    if i > 1e2
        error('error generating the initial size distribution(negative vales: 1st level)!\n')
    end
    i = i + 1;
    if n_pp(1) > 1
        pp_n = round(lognrnd(log(n_pp(1)), log(n_pp(2)), [n_par,1]));
            % Internal primary number distribution
    end
    pp_n(pp_n < 1) = 1; % Minimum attainable number of primaries within...
        % ...each particle
end

% Reinitializing the error generation indices
chk = 1;
i = 0;

% Assigning the internal size distribution to the primaries
while ~ isempty(chk)
    % Generating random lognormal distribution for internal polydispersity
    dist_ip = cell(n_par,1);
    for j = 1 : n_par
        if pp_n(j) > 1
            dist_ip{j} = lognrnd(log(1), log(d_pp(3)), [pp_n(j),1]);
        else
            dist_ip{j} = 1;
        end
    end
    dist_ip = cat(1,dist_ip{:});

    pp_d = repelem(pp_d, pp_n, 1) .* dist_ip;
    chk = find(pp_d <= 0,1);
    if i > 1e2
        error('error generating the initial size distribution (negative vales: 2nd level)!\n')
    end
    i = i + 1;    
end

end

