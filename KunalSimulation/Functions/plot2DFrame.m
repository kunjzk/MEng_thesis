function plot2DFrame(sensor_nH, chip, time_step, frame_number)
    
    frame_to_plot = sensor_nH(:, frame_number);
    n_prot = max(frame_to_plot);
    proton_map = zeros(chip.N_y, chip.N_x);
    time = frame_number*time_step;
    
    for j = 1:chip.N_y
        proton_map(j,:) = frame_to_plot((j-1)*chip.N_x+1:(j-1)*chip.N_x + chip.N_x)';
    end
        
    f3 = figure();
    [X,Y] = meshgrid(1:chip.N_x, 1:chip.N_y);
    h = surf(X, Y, proton_map);

    colorbar;
    title({'Protons on ISFET array, time = ' time 'seconds'}, 'FontSize', 14);
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 chip.N_x]); ylim([1 chip.N_y]); zlim([0 n_prot]);
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;

end

