function plotClustersAndProtons(clusters, protons, chip, solution_height, t)
    figure
    
    
    %Scatter plot DNA clusters
    dna = scatter3(cat(1,clusters.centre_x./chip.isfet_length),...
                   cat(1,clusters.centre_y./chip.isfet_width),...
                   cat(1,clusters.centre_z./chip.height_unit), 'LineWidth', 3);
    
    hold on
    pos_x = (protons(1,:)./chip.isfet_length);
    pos_y = (protons(2,:)./chip.isfet_width);
    pos_z = (protons(3,:)./(chip.height_unit));
    prot = scatter3(pos_x, pos_y, pos_z, '.', 'LineWidth', 3);
    
    
    [X,Y] = meshgrid(1:chip.N_x, 1:chip.N_y);
    h = surf(X, Y, zeros(chip.N_y, chip.N_x), 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
    hold off
    
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 chip.N_x]); ylim([1 chip.N_y]); zlim([0 ceil(solution_height/chip.height_unit)]);
    title({'Release of Protons in 3D Space, time = ' t/60 'minutes'}, 'FontSize', 20);
    % Aesthetics...
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;

end

