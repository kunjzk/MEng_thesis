%Plot initial distribution of DNA clusters
function PlotClustersInitial(num_of_mol, clusters_init, sensor_xsize, sensor_ysize, solution_height, chip)

    figure;
    
    %Scatter plot DNA clusters
    hh = scatter3(cat(1,(clusters_init.centre_x./chip.isfet_length)),...
                    cat(1,(clusters_init.centre_y./chip.isfet_width)),...
                    cat(1,(clusters_init.centre_z./chip.height_unit)), 'LineWidth', 3);
    
    %Add chip layout as RGB image to z=0 plane

    %h = surface([0, 0; sensor_xsize, sensor_xsize], [0, sensor_ysize; 0, sensor_ysize], [0, 0; 0, 0], ones(),'facecolor', 'texturemap', 'edgecolor', 'none');
    %Add labels
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 chip.N_x]); ylim([1 chip.N_y]); zlim([0 ceil(solution_height/chip.height_unit)]);
    title({'Starting Position of Initial Molecules in 3D Space';...
          ['Number of Molecules = ' num2str(num_of_mol)]}, 'FontSize', 20);
    
    % Aesthetics...
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;
    
    hold on
    [X,Y] = meshgrid(1:chip.N_x, 1:chip.N_y);
    h = surf(X, Y, zeros(chip.N_y, chip.N_x), 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
    hold off
    
end