% This script is written to test and monitor the performance and outputs...
    % ...of the program. The post-processing parts, as a result, are not...
    % ...the most efficient in terms of computational cost.

clc
clear
close all

%% Defining the constant physical properties

Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.8e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass';...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};

props = table(Name, Value, Unit, Description); % Creating the table...
    % ...of properties

clear Name Value Unit Description

%% Importing the user-defined parameters

% Extracting data from the input file
file_id = fopen('inputs\MCEM_Params.txt','r');
params = textscan(file_id,'%f%s%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(file_id);

params_val = cell2mat(params(1));
dom_size = params_val(1:3); % Domain dimensions (m)
n_par = params_val(4); % Total number of initial particles (-)
n_pp = params_val(5:6); % Number distibution parameters of primaries (-)
d_pp = params_val(7:9); % Size distribution parameters of primaries (m)
temp_f = params_val(10); % Flow temperature (k)
v_f = params_val(11:13); % Flow velocity vector (m/s)
p_f = params_val(14); % Flow pressure (pa)

% NOTE: See the description tab in the input file for more info on these...
    % ...parameters.

clear file_id params params_val 

%% Initializing the computational domain parameters

% Defining the domain structure (for possible later uses)
dom = struct('size',[]);
dom.size = dom_size;
% Input is the domain size array.

% Loading the fluid properties
fl = struct('temp',temp_f,'v',v_f,'p',p_f,'mu',[],'lambda',[]);
% fl is the fluid information structure for particle-fluid interactions.
[fl.mu, fl.lambda] = CPL.KINETIC(fl,props); % fluid viscosity and...
    % ...mean free path

% Declaring the particle structure
par = struct('pp',[],'n',[],'d',[],'r',[],'v',[],...
'delt',[],'rho',[],'tau',[],'diff',[],'lambda',[],'nnl',[]);
% Inputs are list of primaries characteristics (index, size,...
    % ...and coordinates), number of primary within particles,...
    % ...particles' equivalent position and velocity,...
    % ...their diffusive properties (motion timestep, density,...
    % ...relaxation time, diffusion coefficient, and mean free path),...
    % ...and nearest neighbor list.
% NOTE: Element rows correspond to different aggregates info.

% Calculating the primary particle size and number distributions
[pp_d, par.n] = PAR.INITDIAM(n_par,n_pp,d_pp);

% Initializing the primary particle field; Assigning the indices and sizes
par.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], par.n);

% Generating the initial aggregates (if applicable)
[par.pp, par.d] = PAR.INITMORPH(par.pp, par.n);

% Assigning the particle initial locations
par = PAR.INITLOC(dom_size, par);

% Assigning the particle initial velocities
par.v = PAR.INITVEL(par.pp, par.n, fl.temp, props);

disp("The computational domain is successfully initialized...")

fig_init = figure(1);
UTILS.PLOTPAR(dom_size, par, 1);

fig_nn = figure(2);
UTILS.PLOTNN(dom_size, par, randperm(n_pp,1), 10);

%% Solving equation of motion for the particles

k_max = 50;  % Marching index limit
time = zeros(k_max,1);
fig_anim = figure(3);
t_plt = 1;  % Defining a plotting timeframe criterion

prompt = 'Do you want the animation to be saved? Y/N: ';
str = input(prompt,'s');
if (str == 'Y') || (str == 'y')
    video_par = VideoWriter('outputs\Animation_DLCA.avi'); 
    % Initializing the video file
    video_par.FrameRate = 5;  % Setting frame rate
    open(video_par);  % Opening video file
end

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]);  % Initializing textbar
UTILS.TEXTBAR([1, k_max]);  % Indicating start of marching

for k = 2 : k_max
    
    [par, delt] = MOV.MARCH(par, fl);  % Solving equation of motion
    par.r = MOV.PBC(dom_size, par.r);  % Applying periodic...
        % ...boundary conditions
    
%     par = COL.GROW(par);  % Checking for collisions and updating particle
%     % structures upon new clusterations
    
    if mod(k-2,t_plt) == 0
        UTILS.PLOTPAR(dom_size, par, 0);  % Plotting every t_plt time steps
        drawnow;  % Drawing the plot at the desired time steps
        pause(0.1);  % Slowing down the animation speed
        if (str == 'Y') || (str == 'y')
            frame_now = getframe(fig_anim);  %  Capturing current frame
            writeVideo(video_par, frame_now);  %  Saving the video
        end
    end
    
    time(k) = time(k-1) + delt; % Updating time
    
    UTILS.TEXTBAR([k, k_max]);  % Updating textbar
 
end

if (str == 'Y') || (str == 'y')
    close(video_par);  % Closing the video file
end

