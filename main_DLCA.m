% "main_test" is a script to examine and monitor the performance and...
    %   ...outputs of "MCEM" program.

clc
clear
clf('reset')
close all

%% Part 1. Setting the initial domain conditions

% Setting the run-time record sheet
timetable = struct('start', [], 'end', [], 'total', [],...
    'prerender', [], 'postrender', [], 'render', [],...
    'pregrowth', [], 'postgrowth', [], 'growth', []);
%     'preoverlap', [], 'postoverlap', [], 'overlap', [],...
%     'preconnect', [], 'postconnect', [], 'connect', [],...
%     'preunite', [], 'postunite', [], 'unite', []);
timetable.start = clock;

[params_ud, params_const] = TRANSP.INIT_PARAMS('MCEM_DLCAParams');
    % Initializing the physical parameters to be used in the simulations

[pars, fl] = TRANSP.INIT_DOM(params_ud, params_const); % Initializing...
    % ...the fluid and particle transport variables

% Calculating the primary particle size and number distributions
[pp_d, pars.n] = PAR.INIT_DIAM(params_ud.Value(5), params_ud.Value(6:7),...
    params_ud.Value(8:10));

% Initializing the primary particle field; Assigning the indices and sizes
pars.pp = mat2cell([(1:size(pp_d))', pp_d, zeros(size(pp_d,1),3),...
    (1:size(pp_d))'], pars.n);

% Generating the initial aggregates (if applicable)
pars.pp = PAR.INIT_MORPH_RAND(pars.pp);

% Assigning the particle initial locations
[pars, params_ud] = PAR.INIT_LOC(pars, params_ud);

% Finding the equivalent particle sizes
pars = PAR.SIZING(pars);

% Finding the initial mobility propeties of particles
pars = TRANSP.MOBIL(pars, fl, params_const);

% Assigning the particle initial velocities
pars.v = PAR.INIT_VEL(pars.pp, pars.n, fl.temp, params_const);

% Finding the initial nearest neighbors
ind_trg = (1 : params_ud.Value(5))'; % Indicices of target particles
coef_trg = 5 .* ones(params_ud.Value(5), 1); % Neighboring enlargement...
    % ...coefficients
pars.nnl = COL.NNS(pars, ind_trg, coef_trg);

% Converting the global structure to aggregate objects
% pars = PAR.PAR2AGG(pars);

% Saving the initial population-based particle properties
parsdata = UTILS.SAVEPARS(pars, 0, 1, params_ud);

disp("The computational domain is successfully initialized...")

% Visualizing the initial particle locations, velocities, and nearest...
    % ...neighbor lists
% figure
% h0_pose = UTILS.PLOTPARS(pars, params_ud.Value(1:3),...
%     'equivalent_size', 'on', 'velocity_vector', 'on', 'render', 'on');
% 
% figure
% ind_trg_test = (randperm(params_ud.Value(4),...
%     min([params_ud.Value(4), 2])))'; % A random portion of target...
%         ...indices for nearest neighbor testing
% h0_nntest = UTILS.PLOTNN(pars, params_ud.Value(1:3), ind_trg_test,...
%     coef_trg(ind_trg_test));
% 
% timetable.prerender = clock;
% figure
% UTILS.RENDER(pars);
% timetable.postrender = clock;

%% Part 2: Simulating the particle aggregations

k_max = 1e5; % Marching index limit
time = zeros(k_max,1);
t_rec = 5e2; % Data recording timeframe
% t_plt = 10; % Particle movements plotting ~
% t_nns = 10; % Nearest neighbor search ~
n_agg = zeros(k_max,1); % Real-time number of aggregates
if isa(pars, 'AGG')
    n_agg(1) = length(pars);
else
    n_agg(1) = size(pars.n, 1);
end

% prompt = 'Do you want the aggregation animation to be saved? Y/N: ';
% str = input(prompt,'s'); % Animation saving variable (yes/no)
% str = 'y';
% if (str == 'Y') || (str == 'y')
%    % Checking the existence of animation directory
%     if exist('outputs', 'dir') ~= 7
%         mkdir('outputs');
%     end
%     mkdir('outputs');
%     video_par = VideoWriter('outputs\DLCA_anim.avi'); 
%         % Initializing the video file
%     video_par.FrameRate = 3; % Setting frame rate
%     video_par.Quality = 100; % Setting quality
%     open(video_par); % Opening video file
% elseif (str ~= 'N') && (str ~= 'n')
%     error('Error saving the animation (invalid user response)!\n')
% end

disp('Simulating:');
UTILS.TEXTBAR([0, k_max]); % Initializing textbar

% % Initializing the animation plot
% figure
% h_anim = UTILS.PLOTPARS(pars, params_ud.Value(1:3));
UTILS.TEXTBAR([1, k_max]); % Indicating start of marching

k = 2; % Iteration index
while (k <= k_max) && all(cat(1, pars.n) < (params_ud.Value(5) / 10))
    % Checking if the number of aggregates within the domain is reasonable
    
    [pars, delt] = TRANSP.MARCH(pars, fl, params_const); % Solving...
        % ...equation of motions
    
    pars = TRANSP.PBC(params_ud.Value(2:4), pars); % Applying...
        % ...periodic boundary conditions
    
%     timetable.pregrowth = [timetable.pregrowth; clock];
    pars = COL.GROW(pars); % Checking for...
        % ...collisions and updating particle structures upon new...
        % ...clusterations
%     timetable.postgrowth = [timetable.postgrowth; clock];
    
    pars = PAR.SIZING(pars); % Updating the size-related properties
    
    pars = TRANSP.MOBIL(pars, fl, params_const); % Updating the...
        % ...mobility propeties
    
    %     if mod(k-1, t_nns) == 0
    %         par.nnl = COL.NNS(par, ind_trg, coef_trg); % Finding the...
    %             % ...nearest neighbors
    %     end
    
    % Recording total number of aggregates
    if isa(pars, 'AGG')
        n_agg(k) = length(pars);
    else
        n_agg(k) = size(pars.n, 1);
    end
    
    time(k) = time(k-1) + delt; % Updating time
    
    if mod(k-1, t_rec) == 0
        parsdata = UTILS.SAVEPARS(pars, time, k, params_ud, parsdata);
        % Saving particle data over time
    end
    
    %     if mod(k-1, t_plt) == 0
    %         h_anim = UTILS.PLOTPARS(pars, params_ud.Value(1:3)); % Plotting...
    %             % ...every t_plt time steps
    %         drawnow; % Drawing the plot at the desired time steps
    %         pause(0.1); % Slowing down the animation speed
    %         if (str == 'Y') || (str == 'y')
    %             framenow = getframe(h_anim, [0, 0, 1000, 892.1]);
    %                 % Capturing current frame
    %             writeVideo(video_par, framenow); % Saving the video
    %         end
    %     end
    
    UTILS.TEXTBAR([k, k_max]); % Updating textbar
    
    k = k + 1;
    
end

k = k - 1;

% if (str == 'Y') || (str == 'y')
%     close(video_par); % Closing the video file
% end

%% Part 3: Postprocessing the results

% Morphology of the final population
% timetable.prerender = [timetable.prerender; clock];
% figure
% UTILS.RENDER(pars);
% timetable.postrender = [timetable.postrender; clock];

% Obtaining fractal properties
[df_compiled, kf_compiled] = UTILS.PLOTFRACTALITY(parsdata);

% Plotting kinetic properties
figure
UTILS.PLOTKINETICS(parsdata);

figure
opts_kin.visual = 'on';
if k < k_max
    n_agg(k + 1 : end) = [];
    time(k + 1 : end) = [];
end
[beta, tau, z] = TRANSP.KINETIC(n_agg, time, opts_kin);

% Finalizing the run-time results
% timetable.end = clock;

% timetable.total = UTILS.TIMEDIFF(timetable.start, timetable.end);
% 
% timetable.render = sum(UTILS.TIMEDIFF(timetable.prerender,...
%     timetable.postrender));
% timetable.render = [timetable.render, timetable.render / timetable.total];
% 
% timetable.growth = sum(UTILS.TIMEDIFF(timetable.pregrowth,...
%     timetable.postgrowth));
% timetable.growth = [timetable.growth, timetable.growth / timetable.total];
% 
% timetable.overlap = sum(UTILS.TIMEDIFF(timetable.preoverlap,...
%     timetable.postoverlap));
% timetable.overlap = [timetable.overlap, timetable.overlap /...
%     timetable.growth(1)];
% 
% timetable.connect = sum(UTILS.TIMEDIFF(timetable.preconnect,...
%     timetable.postconnect));
% timetable.connect = [timetable.connect, timetable.connect /...
%     timetable.growth(1)];
% 
% timetable.unite = sum(UTILS.TIMEDIFF(timetable.preunite,...
%     timetable.postunite));
% timetable.unite = [timetable.unite, timetable.unite /...
%     timetable.growth(1)];
