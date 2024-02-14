function [params_ud, params_const] = INIT_PARAMS(fname)
% "INIT_PARAMS" computes various equivalent sizes proposed for...
%   ...fractal aggregates in classic DLCA regime.
% ----------------------------------------------------------------------- %
%
% Output:
%   fname: filename string
%   params_ud: Table of constant parametrs
%   params_const: ~ user-defined parameters
% ----------------------------------------------------------------------- %

% Defining the constant physical parameters used in the soot aggregation...
    % ...problem
Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.86e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass';...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};
params_const = table(Name, Value, Unit, Description); % Creating the...
    % ...table of parameters

% Importing the user-defined parameters of domain, particle and flow

% Extracting user data from the input file
f_adrs = strcat('inputs\', fname, '.txt');
f_id = fopen(f_adrs,'r');
import_data = textscan(f_id,'%s%f%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(f_id);

Name = import_data{1,1};
Value = import_data{1,2};
Unit = import_data{1,3};
Description = import_data{1,4};
params_ud = table(Name, Value, Unit, Description); % The table of...
    % ...user-defined parameters
    % NOTE: See the "Description" column in the input file for more info...
        % ...on these parameters.

end

