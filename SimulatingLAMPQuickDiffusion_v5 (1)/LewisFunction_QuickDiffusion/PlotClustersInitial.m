%Plot initial distribution of DNA clusters
function PlotClustersInitial(num_of_mol, clusters_init, sensor_xsize, sensor_ysize, solution_height, chip_layout)

    figure(7);
    
    %Scatter plot DNA clusters
    if num_of_mol ~= 0
        hh = scatter3(cat(1,clusters_init.centre_x),...
                    cat(1,clusters_init.centre_y),...
                    cat(1,clusters_init.centre_z), 'LineWidth', 3);
    end
    
    %Add chip layout as RGB image to z=0 plane
    hold on
    h = surface([0, 0; sensor_xsize, sensor_xsize], [0, sensor_ysize; 0, sensor_ysize], [0, 0; 0, 0], chip_layout,'facecolor', 'texturemap', 'edgecolor', 'none');
    alpha(h, 0.5);
    hold off

    %Add labels
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([0 sensor_xsize]); ylim([0 sensor_ysize]); zlim([0 solution_height]);
    title({'Starting Position of Initial Molecules in 3D Space';...
          ['Number of Molecules = ' num2str(num_of_mol)]}, 'FontSize', 20);
    
    % Aesthetics...
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;
end