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
pp_n = zeros(n_par,1); % declaring the number array
chk = find(pp_d <= 0,1); % looping criterion
i = 0; % error generation index

% Generating a normally distributed series of random diameters and...
% ...internal primary numbers.
while ~ isempty(chk)
    pp_d = d_pp(2).*randn(n_par,1) + d_pp(1); % Global size distribution
    chk = find(pp_d <= 0,1); % check all the diameters to be positive
    if i > 1e2
        error('error generating the initial size distribution(negative vales: 1st level)!\n')
    end
    i = i + 1;
    pp_n = round(n_pp(2).*randn(n_par,1) + n_pp(1)); % Internal primary...
    % ...number distribution
    pp_n(pp_n < 1) = 1; % Minimum attainable number of primaries within...
    % ...each particle
end

% Random distribution generation for internal polydispersity
dist2_cl = cell(n_par,1);
for j = 1 : n_par
    dist2_cl{j} = randn(pp_n(j),1);
end
dist2 = cat(1,dist2_cl{:});

% Reinitializing the error generation indices
chk = 1;
i = 0;

% Assigning the size distribution to internal aggregate primaries
while ~ isempty(chk)
    pp_d = repelem(pp_d, pp_n, 1) + d_pp(3) .* dist2;
    chk = find(pp_d <= 0,1);
    if i > 1e2
        error('error generating the initial size distribution (negative vales: 2nd level)!\n')
    end
    i = i + 1;    
end

end

