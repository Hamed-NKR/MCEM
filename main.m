clc
clear
close all

%%%%%%


%% Program description

% The present code generates soot fractal aggregates based on...
% diffusion-limited cluster-cluster aggregation (DLCA) approach.

% The numerical algorithm includes transport and aggloneration of...
% nascent primary particles using discrete-element modeling (DEM)...
% in a Langevin dynamics framework.

% Coding by: Hamed Nikookar, PhD candiate in Mechanical Engineering,
% Aerosol lab, University of British Columbia (UBC)

%% Importing the user-defined parameters

% Extracting data from the input file
file_id = fopen('inputs\MCEM_Params.txt','r');
params = textscan(file_id,'%f%s%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(file_id);

params_val = cell2mat(params(1));
dom_size = params_val(1:3); % Domain dimensions (m)
n_pp = params_val(4); % Number of primaries
d_pp = params_val(5:6); % Mean size and standard deviation of primaries (m)
temp_f = params_val(7); % Flow temperature (k)
v_f = params_val(8:10); % Flow velocity vector (m/s)
p_f = params_val(11); % Flow pressure (pa)

%% Initializing the computational domain parameters

% Defining the domain structure
dom = struct('size',[]);
dom.size = dom_size;
% Input is the domain size array.

% Loading the fluid properties
fl = struct('temp',temp_f,'v',v_f,'p',p_f);
% fl is the fluid information structure for particle-fluid interactions.

% Declaring the primary particles structure
pp = struct('ind',[1:n_pp]','ind_agg',zeros(n_pp,1),'d',[],...
    'r',[],'v',[]);
% Inputs are primary particle index, corresponding aggregate index, and...
% position and velocity of the primaries.
% pp rows correspond to different primaries.

% Assigning the primary particle diameters
pp.d = PP.INIT.DIAM(n_pp,d_pp);

% Assigning the primary particle initial locations
pp.r = PP.INIT.LOC(dom_size,n_pp,pp.d);

% Assigning the primary particle initial velocities
pp.v = PP.INIT.VEL(n_pp,pp.d,fl.temp);

disp("The computational domain is successfully initialized...")

fig_pp_init = figure(1);
VIS.STCPLOT(dom.size,pp)

%% Solving equation of motion for the particles

k_max = 10000; % Marching index limit
time = zeros(k_max,1);
rho_pp = 1.8 * 10^3; % Primary particles density ~ Black carbon's bulk density
del_t = (10^-9) * MOV.PROPS(min(pp.d),rho_pp,fl.temp,fl.p); % Marching time step

fig_pp_anim = figure(2);

for k = 2 : k_max
    VIS.DYNPLOT(dom.size,n_pp,pp.d,pp.r)
    pause(0.001)
    for i = 1 : n_pp
        time(k) = time(k-1) + del_t;
        tau_pp = MOV.PROPS(pp.d(i),rho_pp,fl.temp,fl.p);
        [r_pp_new, v_pp_new] = MOV.MARCH(pp.r(i,:),pp.v(i,:),pp.d(i),...
            rho_pp,tau_pp,fl.temp,del_t);
        pp.r(i,:) = r_pp_new;
        pp.v(i,:) = v_pp_new;
    end
end
