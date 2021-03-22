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

% CHANGE: This did not need to be in the loop.
time = 0:del_t:(k_max * del_t);  % times to be looped over

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % initialize textbar
UTILS.TEXTBAR([1, k_max]);  % indicate initial postions are set
for k = 2 : k_max
    
    % CHANGE: Vectorized the particle loop. The program is extremely fast
    % is it doesn't have to plot. 
    
    % TO DO: There are probably too many inputs here. Can this be simplified?
    tau_pp = MOV.PROPS(pp.d,rho_pp,fl.temp,fl.p);  % get particle properties
    [pp.r, pp.v] = MOV.MARCH(pp.r,pp.v,pp.d,...
    	rho_pp,tau_pp,fl.temp,del_t);
    % CHANGE: Assign r and v directly. The multi-step is slower and
    % wastes memory. 
    
    % CHANGE: You should plot after you update the results, or the last
    % step is not plotted. Especially given the VIS.STCPLOT at the
    % beginning plots the original points.
    % CHANGE: Plot every t_plot time steps. Change to t_plot = 1 to plot
    % every step (much slower).
    t_plot = 10;
    if mod(k-1, t_plot)==0
        VIS.DYNPLOT(dom.size,pp.d,pp.r)
        drawnow;  % draw the plot each time step
    end
    
    % CHANGE: I added this textbar utility to indicate progress in command
    % line.
    UTILS.TEXTBAR([k, k_max]);  % update textbar
    
    % TO DO: The particles seem to disappear over time. Likely a bug in
    % applying periodic boundary conditions. 
end
