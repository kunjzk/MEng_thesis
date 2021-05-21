function [ Z ] = sumSensor2(clusters, N)
    % Place a proton on every array for accumarray to output correct dim
    [p,q] = meshgrid(1:N+1, 1:N+1);
    
    x = [clusters.protons_x p(:)'];  %put all the x together since you don't care which cluster they belong to for this function
    y = [clusters.protons_y q(:)'];  %same with y
    Z = accumarray(floor([x', y']), 1) - 1;
    Z = Z(1:N,1:N)';
end