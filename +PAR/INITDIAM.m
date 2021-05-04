function pp_d = INITDIAM(n_pp,d_pp)
% "INITDIAM" assigns a normal size distribution to the primary particles.

% n_pp is the number of primaries.
% d_pp is a 2-member array (2*1) containing mean size and...
% standard deviation of primaries.
% pp_d the diamter size array.

pp_d = zeros(n_pp,1); % declaring the diameter array
chk = find(pp_d <= 0,1); % looping criterion
i = 0; % error generation index

% Generating a normally distributed series of random diameters
while ~ isempty(chk)
    pp_d = d_pp(2).*randn(n_pp,1) + d_pp(1);
    chk = find(pp_d <= 0,1); % check all the diameters to be positive
    if i > 1e3
        error("error generating the initial size distribution!\n")
    end
    i = i + 1;
end

