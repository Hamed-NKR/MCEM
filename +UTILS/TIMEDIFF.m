function delt = TIMEDIFF(t_start, t_end)
% "TIMEDIFF" calculates the difference between two sets of time moments
% ----------------------------------------------------------------------- %
%
% Inputs:
%   t_start: Array of start times
%   t_end: ~ end times
% ----------------------------------------------------------------------- %
%
% Output:
%   delt: Time differences
% ----------------------------------------------------------------------- %

if isempty(t_start) || isempty(t_end)
    delt = 0;
    return
end

m = zeros(12,3); % Initializing month adjustment factor
m(:,1) = (1 : 12)'; % Month indices
% Number of days in each month
m((m(:,1) ==  1) | (m(:,1) == 3) | (m(:,1) == 5) | (m(:,1) == 7) |...
    (m(:,1) == 8) | (m(:,1) == 10) | (m(:,1) == 12), 2) = 31; 
m(m(:,1) ==  2, 2) = 28;
m((m(:,1) ==  4) | (m(:,1) == 6) | (m(:,1) == 9) | (m(:,1) == 11), 2) = 30; 
for i = 1 : 12
    m(i,3) = sum(m(1 : i-1, 2));
end

% Computing the time period between start and end
delt = 31557600 * (t_end(:,1) - t_start(:,1)) +... % Change in years
    86400 * (m(t_end(:,2),3) + t_end(:,3) -...
    m(t_start(:,2),3) - t_start(:,3)) +... % months & days
    3600 * (t_end(:,4) - t_start(:,4)) +... % hours
    60 * (t_end(:,5) - t_start(:,5)) +... % minutes
    t_end(:,6) - t_start(:,6); % seconds

end

