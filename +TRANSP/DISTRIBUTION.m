function [du_dlogd, d_discrete] = DISTRIBUTION(u, d, n)
% "DISTRIBUTION" gives the distribution of a variable with respect to...
%   ...another variable.
% ----------------------------------------------------------------------- %
%
% Inputs:
%   u: Variable to be distributed (N*1 vector)
%   d: Distributor variable (N*1 vector)
%   n: Number of data bins
% ----------------------------------------------------------------------- %
%
% Output:
%   du_dlogd: The (natural) logarithmic distribution of variable "u"...
%       ...with respect to variable "v"
%   d_disctetize: Variable "d" descretized in ten equally sized grid points
% ----------------------------------------------------------------------- %

% Initializing the bin counts
if ~exist('n', 'var'); n = []; end
if isempty(n); n = 10; end

% Discretized form of variable "d"
d_discrete = log(1 : (exp(1) - 1) / n : exp(1))'; % Logarithmic...
    % ...discretization
d_discrete = min(d) + 1.001 * range(d) .* d_discrete; % Identifying the...
    % ...size bins; The 0.001 is just to make sure that the maximum of...
    % ..."d" also falls within the distribution.

% Discretized form of variable "u"
du_dlogd = zeros(n, 1);

for i = 1 : n
    ind = d >= d_discrete(i) & d < d_discrete(i+1); % Indices...
        % ...corresponding to the bins
    du_dlogd(i) = sum(u(ind)) / log(d_discrete(i+1) / d_discrete(i));
        % Assigning derivatives of "u" to the bins
    d_discrete(i) = (d_discrete(i) + d_discrete(i + 1)) / 2;
end
d_discrete(n+1) = [];

end

