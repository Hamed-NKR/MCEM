% "main_test" is a script to examine and monitor the performance and...
    % ...outputs of "MCEM" program.

clc
clear
close all

%% Part 1: Defining the constant physical properties

 % Creating the table of properties
Name = {'rho_bc'; 'M_air'; 'kb'; 'Na'; 'Ru'};
Value = [1.8e3; 28.97e-3; 1.381e-23; 6.022e23; 8.314];
Unit = {'kg/m3'; 'kg/mol'; 'j/k'; 'mol^-1'; 'j/mol.k'};
Description = {'Black Carbon bulk density'; 'Air molar mass';...
    'Boltzmann constant'; 'Avogadro constant'; 'Universal gas constant'};
props = table(Name, Value, Unit, Description);

clear Name Value Unit Description % Cleaning the worksapce

%% Part 2: Importing the user-defined parameters

% Extracting data from the input file
file_id = fopen('inputs\MCEM_Params.txt','r');
import_data = textscan(file_id,'%s%f%s%s','Headerlines',2,...
    'Delimiter','\t','MultipleDelimsAsOne',1);
fclose(file_id);

% Creating the table of parameters
Name = import_data{1,1};
Value = import_data{1,2};
Unit = import_data{1,3};
Description = import_data{1,4};
params = table(Name, Value, Unit, Description);
% NOTE: See the "Description" column for more info on these parameters.

clear file_id import_data ans Name Value Unit Description  % Cleaning...
    % ...the worksapce

%% Part 3: Initializing the computational domain parameters

% Declaring the fluid structure and loading the corresponding imported...
    % ...properties
fl = struct('temp', params.Value(10), 'v', params.Value(11:3),...
    'p', params.Value(14), 'mu', [],'lambda', []);
% "fl" is a structure containing the main physical properties of the...
    % ...background fluid. It is used to consider and study the...
    % ...particle-fluid interactions.
% Fields:
    % temp: Fluid temperature
    % v: ~ velocity
    % p: ~ pressure
    % mu: ~ viscosity
    % lambda: ~ mean free path

[fl.mu, fl.lambda] = CPL.KINETIC(fl,props); % Calculating the fluid...
    % ...viscosity and mean free path

% Declaring the particle structure
par = struct('pp', [], 'n', [], 'd', [], 'r', [], 'v', [], 'm', [],...
'rho', [], 'delt', [], 'tau', [], 'diff', [], 'lambda', [], 'nnl', []);
% "par" is a structure conating the main physical properties of...
    % ...independent particles (whether being single monomers or...
    % ...aggregates) as well as the properties of their constituent...
    % ...primary particles (if applicable).
% Fields:
    % pp: Primary particle characteristics (an N*5 matrix containing...
        % ...their index, size, and location)
    % n: Number of primaries within the particles
    % d: Independent particles equivalent diameter
    % r: ~ spatial location
    % v: ~ velocity
    % m: ~ mass
    % rho: ~ effective density
    % delt: ~ motion time-step
    % tau: ~ relaxation time
    % diff: ~ diffusivity
    % lambda: ~ diffusive mean free path
    % nnl: ~ nearest neighbor list
% NOTE: The rows of each field correspond to characteristics of each...
    % ...independent particle.

% Calculating the primary particle size and number distributions
[pp_d, par.n] = PAR.INITDIAM(params.Value(4), params.Value(5:6),...
    params.Value(7:9));

% Initializing the primary particle field; Assigning the indices and sizes
par.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3)], par.n);

% Generating the initial aggregates (if applicable)
[par.pp, par.d] = PAR.INITMORPH(par.pp, par.n);

% Assigning the particle initial locations
par = PAR.INITLOC(params.Value(1:3), par);

% Assigning the particle initial velocities
par.v = PAR.INITVEL(par.pp, par.n, fl.temp, props);

disp("The computational domain is successfully initialized...")

% Visualizing the initial particle locations and velocities
fig_init = figure(1);
UTILS.PLOTPAR(params.Value(1:3), par,...
     'equivalent_volumetric_size', 'on', 'velocity_vector', 'on');

% Visual testing of the nearest neighbor calculations
fig_nntest = figure(2);
UTILS.PLOTNN(params.Value(1:3), par, randperm(n_pp,1), 10);

%% Part 4: Simulating the particle aggregations

k_max = 50;  % Marching index limit
time = zeros(k_max,1);
fig_anim = figure(3);
t_plt = 1;  % Defining a plotting timeframe criterion

prompt = 'Do you want the aggregation animation to be saved? Y/N: ';
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
    
    [par, delt] = MOV.MARCH(par, fl);  % Solving equation of motions
    par.r = MOV.PBC(params.Value(1:3), par.r);  % Applying periodic...
        % ...boundary conditions
    
%     par = COL.GROW(par);  % Checking for collisions and updating particle
%     % structures upon new clusterations
    
    if mod(k-2,t_plt) == 0
        UTILS.PLOTPAR(params.Value(1:3), par, 0);  % Plotting every...
            % ...t_plt time steps
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

