%% SimulatingLAMPKunal
% Moniri + Lewis with changes
% Castellations not considered
% For simplicity protons can pass through DNA clusters
% No Electrodes
% Trapping regions are used as before

%% Clear workspace and load functions
clear all; 
close all;

addpath('Functions/')
addpath('LewisFunction_QuickDiffusion/')

%% Save Folder

save_D = 'Results/';

%% Flags

%Flag to plot initial distribution of clusters
plotflag_initial = 1;

%Flag for, and frequency of, plotting release of protons
plotflag_proton_release = 1;
plotflag_proton_frequency = 5;

%% Reaction
disp('.')

num_of_mol = 5; % # of initial molecules
solution_height = 1.5E-3; % in metres
speed_of_front = 1.73E-6 ; % metres per second
r_initial = 1.25E-6; % in metres

total_time = 40; % time in minutes
time_step = 5; % in seconds
t = 0:time_step:total_time*60; % in seconds

% Starting time to each cluster is Normal distribution. Values taken from
% fit to real data
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

% Number of protons released is proportional to area of propagation front.
% 'burst' is the proportionality constant. Can use estimateburst script to
% estimate the parameter based on experimental observations.
burst = 5E15;

% Create starting time dist (and sample from it)
pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index 

%Set start time of simulation. signmoid_start(1) to start at the first
%amplification. Enter 1 to start at t=0.
t_start = sigmoid_start(1);

%% Diffusion

%Effective 'rate constant' for H ion diffusion, k_diffusion.
%Value from http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
%for 50 degrees C. Units m^2/s.
diffusivity_H = 1E-8;

%% Chip design
disp('.')

%Design parameters
%Line electrodes have width in x and length in y. Length is equal to the
%sensor_ysize. For no electrodes, best to put electrode_width as zero and halve the
%electrodeisfet_separation
chip.N_x = 78;
chip.N_y = 56; %N_x x N_y array

% isfeet array
chip.isfetisfet_separation = 0E-6; %Separation between the edges of neighbouring ISFETs in the y axis
chip.isfet_width = 20E-6;
chip.isfet_length = 20E-6;
chip.height_unit = 20E-6;
chip.z_scale = 0.75;

%Separation between the edge of the chip and the wall of the reaction
%chamber
chip.wall_separation_xpos = 0e-6;
chip.wall_separation_xneg = 0e-6;
chip.wall_separation_ypos = 0e-6;
chip.wall_separation_yneg = 0e-6;

%% Sensing

% Set the distance above ISFET sensors that capture/detection events will
% occur
detection_distance = 0.1E-6;

%% Trapping
trappingregion_height = 0.5e-3; %all molecules within the box defined by trappingregion_height and trappingregion_chipseparation are trapped
trappingregion_chipseparation = 0.5e-3; %separation between the chip and the trapping region in the x and y directions
trapping_percentage = 50; %trap the the x% of molecules with smallest z coordinate
%% End of Inputs

%% Build Chip

%Calculate the total x and y size of the sensor array
sensor_xsize = chip.N_x*chip.isfet_width + (chip.N_x-1)*(chip.isfetisfet_separation) + chip.wall_separation_xpos + chip.wall_separation_xneg;
sensor_ysize = chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.wall_separation_ypos + chip.wall_separation_yneg;
sensor_vol = sensor_xsize*sensor_ysize*solution_height;

max_cluster_radius = max([sensor_xsize, sensor_ysize, solution_height]);

disp('.')

%% Initialise Clusters

disp('.')
%Clear any existing clusters
if exist('clusters_init','var')
    clear clusters_init
end

%If DNA exists in solution, initialise starting positions of DNA clusters
if num_of_mol ~= 0
    clusters = InitialiseClustersRandom(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial);
end

%Save all variables to file for future reference
param_file = sprintf('%s//Params' , save_D);
if isfile(param_file)
    delete(param_file)
end
save(param_file)

%% Initialise Sensors

%Array to store the number of protons at each sensor over time. col =
%sensor, row = time
sensor_nH = zeros(chip.N_x*chip.N_y, length(t));

%Array to store average nH at sensors over time
sensor_nH_total = zeros(1, length(t));

%% Visualise

%If desired, create a scatter plot to show the distribution of DNA cluster
%centres and the chip layout, as a good reality check
if(plotflag_initial == 1)
    PlotClustersInitial(num_of_mol, clusters, sensor_xsize, sensor_ysize, solution_height, chip)
end

%% Simulation

% Array to store the x,y,z coords of all protons
all_protons = [];
stop_rxn = 0;
tic;
for i = t_start:length(t)-1
    
    if stop_rxn == 1
        disp('stop rxn flag')
        break
    end
    
    protons = [];
    %Carry out next step of LAMP reaction
    num_of_active_sig = sum(i>sigmoid_start);

    for j = 1:num_of_active_sig
        
        % If the radius of any cluster exceeds the maximum permissiable
        % radius (r > len(max dimension)/2), then the cluster is larger
        % than the solution. All reactants have been used up and protons
        % can no longer be produced. This shaves off alot of computation.
        if sum(clusters.radius > max_cluster_radius) > 0
            stop_rxn = 1;
        end
        
        if stop_rxn == 1
            disp('stop rxn flag')
            break
        end
        
        % If cluster is active then update propogation front
        clusters.radius_prev(j) = clusters.radius(j);
        clusters.radius(j) = speed_of_front*(t(i)-t(sigmoid_start(j))); 
        
        % Compute position of randomly distributed protons around propogation front
        numtorelease = round(burst*(volume_of_sphere(clusters.radius(j))  - volume_of_sphere(clusters.radius_prev(j))));

        c = 2*rand(numtorelease,1)-1;
        lon=2*pi*rand(numtorelease,1);
        lat=acos(c);
        a=cos(lon).*sin(lat); %random points on a sphere
        b=sin(lon).*sin(lat);
        
        % These are the produced protons
        additions = [clusters.radius(j)*a'+clusters.centre_x(j);clusters.radius(j)*b'+clusters.centre_y(j); clusters.radius(j)*c'+clusters.centre_z(j)];

        % Remove the computed proton positions outside the array
        additions(:,additions(1,:)>sensor_xsize) = [];
        additions(:,additions(1,:)<0) = [];
        additions(:,additions(2,:)>sensor_ysize) = [];
        additions(:,additions(2,:)<0) = [];
        additions(:,additions(3,:)>solution_height) = [];
        additions(:,additions(3,:)<0) = [];
        
        % Remove protons that fall within other clusters
        for k = 1:num_of_mol
            if(k~=j)
                additions(:,(additions(1,:) - clusters.centre_x(k)).^2 ...
                    + (additions(2,:) - clusters.centre_y(k)).^2 ...
                    + (additions(3,:) - clusters.centre_z(k)).^2 ...
                    <= clusters.radius(k)^2) = [];
            end
        end
        
        protons = [protons additions];
    end
    
    
    
    
    if size(protons, 2) > 0
        
        all_protons = [all_protons protons];
        
        disp('protons released, i=')
        disp(i)
        % Diffusion Equation: Z(t+del_t) = Z(t) + sqrt(2.D.del_t).N(0,1)
        % Rearrange to get del_t, time taken to reach detection height
        del_t = (((detection_distance*ones(1, size(protons,2)) - protons(3,:))./normrnd(0,1, [1, size(protons,2)])).^2)./(2*diffusivity_H);
        
        % Use X, Y diffusion equations to compute final locations of protons
        % protons_final has 4 rows: x,y,z,t
        protons_final = zeros(size(protons, 1), size(protons, 2));
        protons_final(1,:) = protons(1,:) + sqrt(del_t).*(sqrt(2*diffusivity_H)*ones(1, size(protons,2))).*normrnd(0,1, [1, size(protons,2)]);
        protons_final(2,:) = protons(2,:) + sqrt(del_t).*(sqrt(2*diffusivity_H)*ones(1, size(protons,2))).*normrnd(0,1, [1, size(protons,2)]);
        protons_final(3,:) = detection_distance*ones(1, size(protons,2));
        protons_final(4,:) = del_t;
%         disp('protons_final before removal:')
%         disp(protons_final)
    
        % Remove protons outside permissable set of coordinates
        % If a proton hits a wall, remove it. Otherwise we need to consider
        % rebound angles
        protons_final(:,protons_final(1,:)>=sensor_xsize) = [];
        protons_final(:,protons_final(1,:)<=0) = [];
        protons_final(:,protons_final(2,:)>=sensor_ysize) = [];
        protons_final(:,protons_final(2,:)<=0) = [];
        protons_final(:,protons_final(3,:)>=solution_height) = [];
        protons_final(:,protons_final(3,:)<=0) = [];
%         disp('protons_final after removal:')
%         disp(protons_final)
        
        % Figure out which sensor the protons landed on
        % Needs to be 1-indexed, use ceil instead of floow
        protons_final(1,:) = ceil(protons_final(1,:)./chip.isfet_width);
        protons_final(2,:) = ceil(protons_final(2,:)./chip.isfet_length);
        protons_final(4,:) = ceil(protons_final(4,:)./time_step);
        
        % Add protons to that sensor at the right time step        
        for c = 1:size(protons_final, 2)
            %idx = (y-1)*20 + x --> assuming both x and y are 1-indexed
            %because MATLAB
%             x = protons_final(1,c)
%             y = protons_final(2,c)
%             index = x+(y-1)*20
%             timesteps_ = i + protons_final(4,c):length(t)
            sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*20, i+protons_final(4,c):length(t)) = ...
                sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*20, i+protons_final(4,c):length(t)) + 1;
        end
        %% Plot
        
        if plotflag_proton_release ==1
           if mod(i, plotflag_proton_frequency) == 0
              plotClustersAndProtons(clusters, protons, chip, solution_height, t(i));
           end
        end    
        
    else
        disp('zero protons, i=')
        disp(i)
    end
    

end

disp('end of reaction')
disp(toc)

figure
sensor_nH_total = sum(sensor_nH, 1);
plot(t, sensor_nH_total)
title('Total Number of Protons Against Time')

figure
plot(t, sensor_nH)
title('Total Number of Protons per sensor Against Time')

figure
plot(t, mean(sensor_nH,1))
title('Mean Number of Protons Against Time')
