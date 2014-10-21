function [accelerations] = Find_Acceleration(position, parameters)
    
    % calculate the number of planets--the parameters array starts with G and then the masses.
    num_planets = length(parameters) - 1;
    G  = parameters(1);
    mass = parameters(2:10);
    
    % initialize return
    accelerations = zeros(num_planets, 2);
    
    % for each planet 
    for J = 1:num_planets
        % for each other planet
        for K = 1:num_planets
            if J ~= K
                % calculate the distance and direcction
                x_dist = position(K, 1) - position(J, 1);
                y_dist = position(K, 2) - position(J, 2);
                dist_tot = sqrt((x_dist ^ 2) + (y_dist ^ 2));
                x_dir = x_dist / dist_tot;
                y_dir = y_dist / dist_tot;
                % acceleration is the accumulated acceleration in that dir plus the new
                accelerations(J, 1) = accelerations(J, 1) + x_dir * (G * mass(K)/(dist_tot).^2);
                accelerations(J, 2) = accelerations(J, 2) + y_dir * (G * mass(K)/(dist_tot).^2);
            end
        end 
    end
end