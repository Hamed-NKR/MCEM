
% WRITEXYZ  Write particle information to an XYZ file.
%  AUTHOR: Timothy Sipkens, 2021-09-03

function [fc] = WRITEXYZ(pars, n, fn)

if~exist('n', 'var'); n = []; end
if isempty(n); n = 1:length(pars.pp); end

if ~exist('fn', 'var'); fn = []; end

fc = [];  % intialize file contents/output

for ii=1:length(n)
    fc = [fc; ...
        [ii .* ones(size(pars.pp{n(ii)}, 1), 1), ...  % aggregate ID
        pars.pp{n(ii)}(:,2) ./ 2, ...  % radius
        pars.pp{n(ii)}(:,3:5)]];  % XYZ
end

dlmwrite(fn, size(fc,1));
dlmwrite(fn, 'MCEM output', '-append', 'delimiter', '');
if ~isempty(fn)
    dlmwrite(fn, fc, '-append', 'delimiter', ' ')
end

end

