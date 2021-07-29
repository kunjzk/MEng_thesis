function plotSpatialVoltageFrame(final_v, chip, time_step, frame_number)

    disp('.')

    max_v = max(max(final_v));
    min_v = min(min(final_v));
    
    frame_to_plot = final_v(:, frame_number);
    n_prot = max(frame_to_plot);
    v_map = zeros(chip.N_y, chip.N_x);
    time = frame_number*time_step;
    
    for j = 1:chip.N_y
        v_map(j,:) = frame_to_plot((j-1)*chip.N_x+1:(j-1)*chip.N_x + chip.N_x)';
    end
        
    f3 = figure();
    [X,Y] = meshgrid(1:chip.N_x, 1:chip.N_y);
    h = surf(X, Y, v_map);

    colorbar;
    title({'Spatial Voltage Output, time = ' time 'seconds'}, 'FontSize', 16);
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 chip.N_x]); ylim([1 chip.N_y]); zlim([min_v max_v]);
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;
    
    
end

