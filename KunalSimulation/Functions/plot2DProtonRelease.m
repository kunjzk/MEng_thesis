function plot2DProtonRelease(sensor_nH, chip, time_step, filename)

    max_protons = max(max(sensor_nH));
    proton_map = zeros(chip.N_y, chip.N_x, size(sensor_nH,2));
    
    % Do it manually because how does reshape work?!
    for i = 1:size(sensor_nH, 2)
        for j = 1:chip.N_y
            proton_map(j,:,i) = sensor_nH((j-1)*chip.N_x+1:(j-1)*chip.N_x + chip.N_x, i)';
        end
        %proton_map(:,:,i) = reshape(sensor_nH(:,i), chip.N_y, chip.N_x);
    end
        
    f3 = figure();
    [X,Y] = meshgrid(1:chip.N_x, 1:chip.N_y);
    h = surf(X, Y, zeros(chip.N_y, chip.N_x));

    colorbar;
    title({'Release of Protons in 3D Space, time = 0 minutes'}, 'FontSize', 20);
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim([1 chip.N_x]); ylim([1 chip.N_y]); zlim([0 max_protons]);
    set(gca,'FontSize',16);
    box('on'); grid on; grid minor;
    
    for i = 1:size(proton_map, 3)
        h.ZData = proton_map(:,:,i);
        title({'Release of Protons in 3D Space, time =' i*time_step/60 'minutes'}, 'FontSize', 20);
        %drawnow limitrate;
        pause(0.1);
        mymov(i) = getframe(f3);
    end
    
    
v = VideoWriter(filename);
open(v);
writeVideo(v, mymov);
close(v);
    

end

