%% SimulatingLAMPKunal
% Moniri + Lewis with changes
% Castellations not considered
% For simplicity protons can pass through DNA clusters
% No Electrodes
% Trapping regions are used as before

%% Clear workspace and load functions
clear all; 
close all;

addpath('/Users/kunal/Desktop/thesis/KunalSimulation/Functions/')

time_steps = [1 3 5 10];

%% Flags

%Flag to plot initial distribution of clusters
plotflag_initial = 1;

%Flag for, and frequency of, plotting release of protons
plotflag_proton_release = 0;
plotflag_proton_frequency = 10;

plotflag_2d_surface = 0;
plot_2d_surface_frequency = 5;
save_2d_surface = 0;

stop_rxn_tolerance = 10;

%% Reaction
disp('.')

num_of_mol = 100; % # of initial molecules
speed_of_front = 1.73E-6 ; % metres per second
r_initial = 1.25E-6; % in metres

% Starting time to each cluster is Normal distribution. Values taken from
% fit to real data
mu_start_time = 13; % in minutes
sigma_start_time = 3; % in minutes

removed_1 = zeros(1, 2401);
removed_3 = zeros(1, 801);
removed_5 = zeros(1, 481);
removed_10 = zeros(1,241);

not_removed_1 = zeros(1, 2401);
not_removed_3 = zeros(1, 801);
not_removed_5 = zeros(1, 481);
not_removed_10 = zeros(1,241);

released_1 = zeros(1, 2401);
released_3 = zeros(1, 801);
released_5 = zeros(1, 481);
released_10 = zeros(1, 241);

saturation_1 = 0;
saturation_3 = 0;
saturation_5 = 0;
saturation_10 = 0;


for counter=1:4
    
    total_time = 40; % time in minutes
    time_step = time_steps(counter); % in seconds
    t = 0:time_step:total_time*60; % in seconds

    % Number of protons released is proportional to area of propagation front.
    % 'burst' is the proportionality constant. Can use estimateburst script to
    % estimate the parameter based on experimental observations.
    burst = 3E13; %1.15E9 for 6mV, 4E8 for pH

    % Create starting time dist (and sample from it)
    pd = makedist('Normal', 'mu', mu_start_time, 'sigma', sigma_start_time);
    sigmoid_start = round(60 * sort(random(pd, 1, num_of_mol)) / time_step); % notice conversion to index 


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
    chip.N_x = 56;
    chip.N_y = 78; %N_x x N_y array

    % isfet array
    % Lets assume isfet separations are 1/5 of isfet dimension
    chip.isfetisfet_separation = 2E-6;
    % Distance between isfets and sensor (both dimensions)
    chip.isfet_startSeparation = 2E-6;
    chip.isfet_endSeparation = 2E-6;
    % ISFET dimensions
    chip.isfet_width = 35E-6;
    chip.isfet_length = 37E-6;
    % For plotting
    chip.height_unit = 20E-6;

    %Calculate the total x and y size of the sensor array
    sensor_xsize = chip.N_x*chip.isfet_width + (chip.N_x-1)*(chip.isfetisfet_separation) + chip.isfet_startSeparation + chip.isfet_endSeparation;
    sensor_ysize = chip.N_y*(chip.isfet_length) + (chip.N_y-1)*(chip.isfetisfet_separation) + chip.isfet_startSeparation + chip.isfet_endSeparation;

    %%

    %Separation between the edge of the sensor array and the wall of the reaction
    %chamber
    chip.wall_separation_xpos = 2E-6;
    chip.wall_separation_xneg = 2E-6;
    chip.wall_separation_ypos = 2E-6;
    chip.wall_separation_yneg = 2E-6;

    reaction_x_size = sensor_xsize + chip.wall_separation_xpos + chip.wall_separation_xneg;
    reaction_y_size = sensor_ysize + chip.wall_separation_ypos + chip.wall_separation_yneg;

    solution_vol = 12*10^(-9); %12uL
    solution_height = solution_vol/(sensor_xsize*sensor_ysize);
    sensor_vol = sensor_xsize*sensor_ysize*solution_height;

    % 1m^3 = 1000L
    avogadro = 6.02E23;
    initial_pH = 8.75;
    initial_conc = 10^(-1*initial_pH); %moles/L
    initial_n_prot = 10^(-1*initial_pH)*avogadro*solution_vol/1000; %molecules
    concentration = ((initial_n_prot/solution_vol)/avogadro)*1000; %moles/L
    initial_pH_reverse = log10(concentration)*-1;

    %% Detection Distance = Debye Length

    debye.perm_space = 8.8544E-12;
    debye.perm_water = 65;
    debye.k = 1.3806E-23;
    debye.T = 336.15;
    debye.N_a = 6.023E23;
    debye.e = 1.6022E-19;
    debye.I = 2;

    % Set the distance above ISFET sensors that capture/detection events will
    % occur
    debye_length = sqrt(debye.perm_space*debye.perm_water*debye.k*debye.T/(2*debye.N_a*debye.e^2*debye.I));
    detection_distance = debye_length;

    %% Permissable proton coordinates

    % Build set of allowable coordinates for protons to be sensed
    chip.x_allowed_start = zeros(1, chip.N_x);
    chip.x_allowed_end = zeros(1, chip.N_x);
    chip.x_allowed_start(1) = chip.wall_separation_xneg + chip.isfet_startSeparation;
    for c = 2:chip.N_x
        chip.x_allowed_start(c) = round(chip.x_allowed_start(c-1) + chip.isfet_width + chip.isfetisfet_separation, 10);
    end
    chip.x_allowed_end = chip.x_allowed_start + chip.isfet_width;

    chip.y_allowed_start = zeros(1, chip.N_y);
    chip.y_allowed_end = zeros(1, chip.N_y);
    chip.y_allowed_start(1) = chip.wall_separation_yneg + chip.isfet_startSeparation;
    for c = 2:chip.N_y
        chip.y_allowed_start(c) = chip.y_allowed_start(c-1) + chip.isfet_length + chip.isfetisfet_separation;
    end
    chip.y_allowed_end = chip.y_allowed_start + chip.isfet_length;

    %% End of Inputs

    %% Build Chip

    max_cluster_radius = sqrt(reaction_x_size^2 + reaction_y_size^2 + solution_height^2);

    %% Initialise Clusters

    %If DNA exists in solution, initialise starting positions of DNA clusters
    if num_of_mol ~= 0
        clusters = InitialiseClustersRandom(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial);
        %clusters = InitializeClusterHeight(num_of_mol, sensor_xsize, sensor_ysize, solution_height, r_initial, 0.9);
    end

    %Set start time of simulation. signmoid_start(1) to start at the first
    %amplification. Enter 1 to start at t=0.
    t_start = sigmoid_start(1);

    %% Initialise Sensors
    %Array to store the number of protons at each sensor over time. col =
    %sensor, row = time
    sensor_nH = zeros(chip.N_x*chip.N_y, length(t));
    cum_nH = zeros(chip.N_x*chip.N_y, length(t));
    total_protons_in_solution = zeros(1, length(t));

    %% Simulation
    % Array to store the x,y,z coords of all protons
    stop_rxn = 0;
    tic;
    steps_with_no_protons_released = 0;
    


    for i = t_start:length(t)-1
        
        total_prot_removed_in_this_step = 0;
        total_prot_free = 0;
        total_prot_released_in_this_step = 0;

        if steps_with_no_protons_released == stop_rxn_tolerance
            stop_rxn = 1;
        end

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
            %numtorelease = round(burst*(clusters.radius(j)^2));

            c = 2*rand(numtorelease,1)-1;
            lon=2*pi*rand(numtorelease,1);
            lat=acos(c);
            a=cos(lon).*sin(lat); %random points on a sphere
            b=sin(lon).*sin(lat);

            % These are the produced protons
            additions = [clusters.radius(j)*a'+clusters.centre_x(j);clusters.radius(j)*b'+clusters.centre_y(j); clusters.radius(j)*c'+clusters.centre_z(j)];
            total_released = size(additions,2);
            % Remove the computed proton positions outside the sensor
            additions(:,additions(1,:)>sensor_xsize+chip.wall_separation_xneg) = []; %this should be xneg
            additions(:,additions(1,:)<chip.wall_separation_xneg) = [];
            additions(:,additions(2,:)>sensor_ysize+chip.wall_separation_yneg) = [];
            additions(:,additions(2,:)<chip.wall_separation_yneg) = [];
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
            
            total_remaining = size(additions,2);

            protons = [protons additions];
            total_prot_removed_in_this_step = total_prot_removed_in_this_step + total_released-total_remaining;
            total_prot_free = total_prot_free + total_remaining;
            total_prot_released_in_this_step = total_prot_released_in_this_step + total_released;
        end

        if size(protons, 2) > 0

            steps_with_no_protons_released = 0;

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

            protons_final = checkCoordinatesAndTime(protons_final, chip, time_step, i, length(t));

            % Figure out which sensor the protons landed on
            % Needs to be 1-indexed, use ceil instead of floow
            protons_final(1,:) = ceil((protons_final(1,:)-(chip.isfet_startSeparation+chip.wall_separation_xneg))./(chip.isfet_width+chip.isfetisfet_separation));
            protons_final(2,:) = ceil((protons_final(2,:)-(chip.isfet_startSeparation+chip.wall_separation_yneg))./(chip.isfet_length+chip.isfetisfet_separation));
            protons_final(4,:) = ceil(protons_final(4,:)./time_step);

            total_protons_in_solution(i) = total_protons_in_solution(i-1) + size(protons, 2)-size(protons_final,2);


            % Add protons to that sensor at the right time step     
    %             x = protons_final(1,c)
    %             y = protons_final(2,c)
    %             index = x+(y-1)*chip.N_x
    %             timesteps = i + protons_final(4,c):length(t) (protons persist
    %             from the time they were released till the end of the
    %             reaction)
            for c = 1:size(protons_final, 2)
                % Instantaneous
                sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c)) = ...
                    sensor_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c)) + 1;

                % Cumulative
                cum_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c):length(t)) = ...
                    cum_nH(protons_final(1,c)+(protons_final(2,c)-1)*chip.N_x, i+protons_final(4,c):length(t)) + 1;
            end

        else
            disp('zero protons, i=')
            disp(i)
            total_protons_in_solution(i) = total_protons_in_solution(i-1);
            %If (1) there is a positive number of protons in the solution 
            % and (2) all clusters have started amplifying,
            % AND 
            % No new protons were released, this is a signal to stop rxn
            if max(total_protons_in_solution) > 0 %&& min(clusters.radius)>r_initial
               steps_with_no_protons_released = steps_with_no_protons_released + 1;
            end
        end
        
        %removed_prot_per_step = [removed_prot_per_step total_prot_removed_in_this_step];
        
        if counter == 1
            removed_1(i) = total_prot_removed_in_this_step;
            not_removed_1(i) = total_prot_free;
            released_1(i) = total_prot_released_in_this_step;
        end
        
        if counter == 2
            removed_3(i) = total_prot_removed_in_this_step;
            not_removed_3(i) = total_prot_free;
            released_3(i) = total_prot_released_in_this_step;
        end
        
        if counter == 3
            removed_5(i) = total_prot_removed_in_this_step;
            not_removed_5(i) = total_prot_free;
            released_5(i) = total_prot_released_in_this_step;
        end
        
        if counter == 4
            removed_10(i) = total_prot_removed_in_this_step;
            not_removed_10(i) = total_prot_free;
            released_10(i) = total_prot_released_in_this_step;
        end
        
    end

    disp('end of reaction')
    disp(toc)
    
    if counter == 1
        t1 = t;
        cum_nH_1 = cum_nH;
        saturation_1 = i;
        %removed_1 = removed_prot_per_step;
    end
    
    if counter == 2
        t3 = t;
        cum_nH_3 = cum_nH;
        saturation_3 = i;
        %removed_3= removed_prot_per_step;
    end
    
    if counter == 3
        t5 = t;
        cum_nH_5 = cum_nH;
        saturation_5 = i;
        %removed_5 = removed_prot_per_step;
    end
    
    if counter == 4
        t10 = t;
        cum_nH_10 = cum_nH;
        saturation_10 = i;
        %removed_10 = removed_prot_per_step;
    end
    
end




%%
frame_rates = [1 3 5 10]

cum_nH_1 = zeros(chip.N_x*chip.N_y, 2400);
cum_nH_3 = zeros(chip.N_x*chip.N_y, 800);
cum_nH_5 = zeros(chip.N_x*chip.N_y, 480);
cum_nH_10 = zeros(chip.N_x*chip.N_y, 240);

for counter = 1:4
    frame_rate = frame_rates(counter);
    frame_times = 0:frame_rate:total_time*60;
    frames_to_save = frame_times + 1;

    if counter == 1
        cum_nH_1 = cum_nH(:, [frames_to_save]);
        t1 = frame_times;
    end

    if counter == 2
        cum_nH_3 = cum_nH(:, [frames_to_save]);
        t3 = frame_times;
    end

    if counter == 3
        cum_nH_5 = cum_nH(:, [frames_to_save]);
        t5 = frame_times;
    end
    
    if counter == 4
        cum_nH_10 = cum_nH(:, [frames_to_save]);
        t10 = frame_times;
    end

end

%% Plot temporal to figure out region of interest
figure
plot(t1, sum(cum_nH_1), 'DisplayName', 'Time step = 1', 'LineWidth', 2)
hold on
plot(t3, sum(cum_nH_3), 'DisplayName', 'Time step = 3', 'LineWidth', 2)
plot(t5, sum(cum_nH_5), 'DisplayName', 'Time step = 5', 'LineWidth', 2)
plot(t10, sum(cum_nH_10), 'DisplayName', 'Time step = 10', 'LineWidth', 2)
hold off
grid on 
grid minor
legend
title('Temporal Signal vs Time, Varying Time Step')
xlabel('Time (s)')
ylabel('Number of Protons Detected by ISFET Array')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(t1, cumsum(removed_1), 'DisplayName', '{\tau} = 1', 'LineWidth', 1)
hold on
plot(t3, cumsum(removed_3), 'DisplayName', '{\tau} = 3', 'LineWidth', 1)
plot(t5, cumsum(removed_5), 'DisplayName', '{\tau} = 5', 'LineWidth', 1)
plot(t10, cumsum(removed_10), 'DisplayName', '{\tau} = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Cumulative Number Protons Removed')
xlabel('Time (s)')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(t1, (not_removed_1), 'DisplayName', '{\tau} = 1', 'LineWidth', 1)
hold on
plot(t3, (not_removed_3), 'DisplayName', '{\tau} = 3', 'LineWidth', 1)
plot(t5, (not_removed_5), 'DisplayName', '{\tau} = 5', 'LineWidth', 1)
plot(t10, (not_removed_10), 'DisplayName', '{\tau} = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Number of Protons Remaining vs Time')
xlabel('Time (s)')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

figure
plot(t1, cumsum(not_removed_1), 'DisplayName', '{\tau} = 1', 'LineWidth', 1)
hold on
plot(t3, cumsum(not_removed_3), 'DisplayName', '{\tau} = 3', 'LineWidth', 1)
plot(t5, cumsum(not_removed_5), 'DisplayName', '{\tau} = 5', 'LineWidth', 1)
plot(t10, cumsum(not_removed_10), 'DisplayName', '{\tau} = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Cumulative Number of Protons Remaining')
xlabel('Time (s)')
ylabel('Number of protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)


%% 

figure
plot(t1, cumsum(released_1), 'DisplayName', '{\tau} = 1', 'LineWidth', 1)
hold on
plot(t3, cumsum(released_3), 'DisplayName', '{\tau} = 3', 'LineWidth', 1)
plot(t5, cumsum(released_5), 'DisplayName', '{\tau} = 5', 'LineWidth', 1)
plot(t10, cumsum(released_10), 'DisplayName', '{\tau} = 10', 'LineWidth', 1)
hold off
grid on 
grid minor
legend
title('Number of New Protons Released vs Time')
xlabel('Time (s)')
ylabel('Number of Protons')
grid on
grid minor
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14)

%%

mean(removed_1)
mean(removed_3)
mean(removed_5)
mean(removed_10)

%% Conversion to voltage

% Constants
threshold.temperature = 336.15; %In Kelvin. 336.15K for 63C for LAMP reaction
threshold.k_B = 1.4E-23; %Boltzmann's constant
threshold.q = -1.602E-19; %Electron charge
threshold.permittivity_water = 65; %Roughly 65 at 63C
threshold.permittivity_vacuum = 8.854E-12;
threshold.d_OHP = 1E-7; %Rough estimate
threshold.n_conc = 1.1*1000*6.02E23; %In mol/L (M). Concentration of charged species
threshold.Ka = 10^(11.9); % mol/m2
threshold.Kb = 10^(-1.9); % mol/m2
threshold.Ns = 5E18; %sites/m^2
avogadro = 6.02E23;

%%

for counter = 1:4
    frame_rate = frame_rates(counter);

    if counter == 1
        cum_nH = cum_nH_1;
        t = t1;
    end

    if counter == 2
        cum_nH = cum_nH_3;
        t = t3;
    end

    if counter == 3
        cum_nH = cum_nH_5;
        t = t5;
    end
    
    if counter == 4
        cum_nH = cum_nH_10;
        t = t10;
    end
    
    
    %% New idea: Corrected ahmad equations, but q=charge is negative

    initial_pH = 8.75;
    initial_psi_0 = -150*10^(-3); %from fig 5 
    initial_bulk_h = 10^(-1*initial_pH)*1000; %mol/m^3
    initial_surface_h = initial_bulk_h*exp(-threshold.q*initial_psi_0/(threshold.k_B*threshold.temperature));
    n_prot = 8; %initial_surface_h*(chip.isfet_length*chip.isfet_width)*avogadro;
    cum_nH_include = cum_nH_1+n_prot;
    isfet_debye_length = 6.57*10^-9; %8.9*10^(-15); %computed debye length, defauly 10^-9

    % Compute pH in mol/L
    conc_H = (cum_nH_include./avogadro)./(chip.isfet_width*chip.isfet_length*isfet_debye_length*1000);
    pH = log10(conc_H).*-1;

    %del_pH = diff(pH,1,2);

    del_v = pH.*(2.3*59*10^(-3)*threshold.k_B*threshold.temperature/threshold.q);

    figure
    plot(t1, del_v)
    hold on
    plot(t1, mean(del_v), 'LineWidth', 3, 'Color', 'Black')
    hold off
    title('Method 9 - correct eqn & q is -ve')


    %% Add drift

    % Tau is the time at which the population is reduced to 0.367 the initial
    % value. From experiments this seems to be ~1000s
    % but I reduced it to make the decay steeper so that the kink in the temporal signal from the
    % positive experiment looks subtler and more realistic 
    mean_tau = 750;
    std_tau = 50;
    tau = makedist('Normal', 'mu', mean_tau, 'sigma', std_tau);

    % Take this from the drift modelling thesis
    mean_beta = 0.722;
    std_beta = 0.050;
    beta = makedist('Normal', 'mu', mean_beta, 'sigma', std_beta);

    % The final value of the output
    % Seems to decay by ~100mV
    mean_Q = -100E-3;
    std_Q = 5E-3;
    Q = makedist('Normal', 'mu', mean_Q, 'sigma', std_Q);

    % Create a unique drift for each sensor
    drift = zeros(chip.N_x*chip.N_y, length(t));
    for i=1:size(drift, 1)
        exp_term = 1-exp(((t./random(tau)).^random(beta).*-1));
        drift(i, :) = exp_term.*random(Q);
    end

    %% Add AWGN

    awgn = makedist('Normal', 'mu', 0, 'sigma', 1);
    awgn_scale_factor = 0.568E-3; % magnitude of AWGN

    for i=1:size(drift, 1)
        whitenoise = random(awgn,[1,length(t)]).*awgn_scale_factor;
        drift(i,:) = drift(i, :) + whitenoise;
    end

    %% Negative Reaction
    cum_nH_neg = zeros(chip.N_x*chip.N_y, length(t));

    initial_pH = 8.75;
    initial_psi_0 = -150*10^(-3); %from fig 5 
    initial_bulk_h = 10^(-1*initial_pH)*1000; %mol/m^3
    initial_surface_h = initial_bulk_h*exp(-threshold.q*initial_psi_0/(threshold.k_B*threshold.temperature));
    n_prot = initial_surface_h*(chip.isfet_length*chip.isfet_width)*avogadro;
    cum_nH_include_neg = cum_nH_neg+n_prot;
    isfet_debye_length = 10^-9; %8.9*10^(-15); %computed debye length, defauly 10^-9

    % Compute pH in mol/L
    conc_H_neg = (cum_nH_include_neg./avogadro)./(chip.isfet_width*chip.isfet_length*isfet_debye_length*1000);
    pH_neg = log10(conc_H_neg).*-1;

    %del_pH = diff(pH,1,2);

    del_v_neg = pH_neg.*(2.3*59*10^(-3)*threshold.k_B*threshold.temperature/threshold.q);


    %% Final plot - TEMPORAL

    final_out = drift + del_v;
    final_neg = drift + del_v_neg;
    % subplot(1,2,1)

%     figure
%     sensor_signal = plot(t, final_out);
%     alpha(.1);
%     hold on
%     m = plot(t, mean(final_out), 'LineWidth', 4, 'Color', 'Black', 'DisplayName', 'Mean');
%     hold off
%     title({'ISFET Array Output (Positive Sample)' 'Frame Rate = ' frame_rates(counter)},'FontSize', 16);
%     xlabel('Time (s)')
%     ylabel('Voltage output (V)')
%     legend(m);
%     ax = gca;
%     ax.YAxis.Exponent = -3;
%     grid on
%     grid minor
%     a = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',a,'fontsize',14)
    
    %% Final output - SPATIAL std
%     figure
%     plot(t, std(final_out))
%     title({'Standard Deviation of final output, frame rate = ' frame_rates(counter)})
%     
%     figure
%     plot(t, std(final_out) - std(final_neg))
%     title({'Diff in std, frame rate = ' frame_rates(counter)})
%     
%     plotSpatialVoltageFrame(final_out, chip, frame_rates(counter), ceil(510/frame_rates(counter)))
%     plotSpatialVoltageFrame(final_out, chip, frame_rates(counter), ceil(511/frame_rates(counter)))

    diffv = [diff(mean(final_out))];
    figure
    plot(1:1:length(diffv), diffv.*-1)
    title({'Difference in average voltage across frames, frame rate = ' frame_rates(counter)})
    
end




%% Final plot - SPATIAL

plotSpatialVoltage(final_out, chip, time_step)


%% Final plot - spatial diff
delv = diff(final_out,1,2);
frame_diff = [zeros(chip.N_x*chip.N_y, 1) delv];
plotSpatialVoltage(frame_diff, chip, time_step)

%% Final plot - cum_nH

plotSpatialVoltage(cum_nH, chip, time_step)

%%
plot2DFrame(cum_nH_1, chip, 1, 700)
plotSpatialVoltageFrame(final_out, chip, 1, 700)
% plot2DFrame(cum_nH_3, chip, 3, 134)
% plot2DFrame(cum_nH_5, chip, 5, 80)

%%
final_neg = drift + del_v_neg;

figure
plot(t1, std(final_neg))
title('Negative output std')

figure
plot(t1, std(final_out))
title('Positive output std')

figure
plot(t1, std(final_out)-std(final_neg))
title('diff in std')

figure
plot(t1, std(cum_nH_1))
title('Proton std')

figure
plot(t1, mean(final_out)-mean(drift))
title('temporal voltage signal')
