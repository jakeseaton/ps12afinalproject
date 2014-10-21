% This script uses the Euler algorithm to trach N_steps of 8 planets around the sun. 
% It utilizes nested for loops to conduct kinematic analysis at each stage.
% I - The current step (out of N). this is the outermost for loop
% J - The current planet
% K - The current other planet
% When we calculate the acceleration of a planet, we iterate through all
% of the other planets and calculate the acceleration toward them, 
% accumulating acceleration in the X and Y direction as we go. 
% There is a bug--When iterating through the other planets we go over the current
% planet, the distance will be zero which will create a division by zero error
% The solution to this is probably a if J = K block around the update
% TO DO: Intial values. Fix divide by zero bug
% Did I forget the squared R?
% Function to compute the acceleartion
% Function that does the midpoint algorithm
% Pass the acceleration function as an input 
% Matlab has built in differential equation solvers.
% ODE 113

clear all

% Steps to iterate
N_steps = 100000;

% physical constants;
G = 6.67E-11;

% Time array
T = 3.15569e7; % one year in seconds
t = linspace(0,2*T, N_steps);
dt = mean(diff(t));

% Initialize arrays with structure array(Planet, Dir, N)
position = zeros(8,2,N_steps);
velocity = zeros(8,2,N_steps);
acceleration = zeros(8,2,N_steps);

% Initialze the Sun
xsun = 0; 
ysun = 0;
msun = 1.989*10^30;
vsun = 0;

% Initialize planets, each with an initial x position, y position, and mass

% Mercury
x01 = 57*10^9;
y01 = 0;
vx01 = 0;
vy01 = 47.9*10^3;
m1 = 328.5*10^21;

% Venus
x02 = 108*10^9;
y02 = 0;
vx02 = 0;
vy02 = 35*10^3;
m2 = 4.867*10^24;

% Earth
x03 = 150*10^9;
y03 = 0;
vx03 = 0;
vy03 = 29.8*10^3;
m3 = 5.972*10^24;

% Mars 
x04 = 228*10^9;
y04 = 0;
vx04 = 0;
vy04 = 24.1*10^3;
m4 = 629*10^21;

% Jupiter
x05 = 779*10^9;
y05 = 0;
vx05 = 0;
vy05 = 13.1*10^3;
m5 = 1.898*10^27;

% Saturn 
x06 = 1.43*10^12;
y06 = 0;
vx06 = 0;
vy06 = 9.6*10^3;
m6 = 568.3*10^24;

% Uranus 
x07 = 2.88*10^12;
y07 = 0;
vx07 = 0;
vy07 = 6.8*10^3;
m7 = 86.8*10^24;

% Neptune
x08 = 4.5*10^12;
y08 = 0;
vx08 = 0;
vy08 = 5.4*10^3;
m8 = 102.2*10^24;

% Mass array
mass = [m1, m2, m3, m4, m5, m6, m7, m8];
x0_array = [x01, x02, x03, x04, x05, x06, x07, x08];
y0_array = [y01, y02, y03, y04, y05, y06, y07, y08];
vx0_array = [vx01, vx02, vx03, vx04, vx05, vx06, vx07, vx08];
vy0_array = [vy01, vy02, vy03, vy04, vy05, vy06, vy07, vy08];

% initialize positions and velocities
for J = 1:8
	position(J, 1, 1) = x0_array(J);
	position(J, 2, 1) = y0_array(J);

	velocity(J, 1, 1) = vx0_array(J);
	velocity(J, 2, 1) = vy0_array(J);
end

% initialize accelerations
% iterate over planets
for J = 1:8


    % the initial acceleration is initialized as the acceleration toward the sun
		x_dist_sun = xsun-position(J,1,1);
		y_dist_sun = ysun-position(J,2,1);
		dist_to_sun = sqrt(x_dist_sun ^ 2 + y_dist_sun ^ 2);
		x_dir_sun = x_dist_sun / dist_to_sun;
		y_dir_sun = y_dist_sun / dist_to_sun;
		acceleration (J, 1, 1) = x_dir_sun * G * msun / ((dist_to_sun).^2);
		acceleration (J, 2, 1) = y_dir_sun * G * msun / ((dist_to_sun).^2);
        
        for K = 1:8
            if J ~= K
            % the new acceleration in the x dir is equal to the previous one plus the new one
            x_dist = position(K, 1, 1) - position(J, 1, 1);
            y_dist = position(K, 2, 1) - position(J, 2, 1);
            dist_tot = sqrt((x_dist ^ 2) + (y_dist ^ 2));
            x_dir = x_dist / dist_tot;
            y_dir = y_dist / dist_tot;
            % they accumulate acceleration each time
            acceleration(J, 1, 1) = acceleration(J, 1, 1) + x_dir * (G * mass(K)/(dist_tot).^2);
            acceleration(J, 2, 1) = acceleration(J, 2, 1) + y_dir * (G * mass(K)/(dist_tot).^2);
            end
        end 
end

% Iterate through N_steps
for I = 2: N_steps	

    % initialize half step arrays
    position_half = zeros(8, 2);
    velocity_half = zeros(8,2);
    acceleration_half = zeros(8,2);

	% position and velocity half step 
	for J = 1:8

        %% position half step
        position_half(J,1) = position(J,1,I-1) + velocity(J,1,I-1)*dt/2;
        position_half(J,2) = position(J,2,I-1) + velocity(J,2,I-1)*dt/2;
		
        % velocity half step
        velocity(J,1) = velocity(J,1,I-1) + acceleration(J,1,I-1)*dt;
        velocity(J,2) = velocity(J,2,I-1) + acceleration(J,2,I-1)*dt;

	end 
    
    % acceleration half step
    for J = 1:8
        % the new acceleration is initialized as the acceleration toward the sun
        x_dist_sun = xsun - position_half(J,1);
        y_dist_sun = ysun - position_half(J,2);
        dist_to_sun = sqrt(x_dist_sun ^ 2 + y_dist_sun ^ 2);
        x_dir_sun = x_dist_sun / dist_to_sun;
        y_dir_sun = y_dist_sun / dist_to_sun;
        acceleration_half(J, 1) = x_dir_sun * G * msun / ((dist_to_sun).^2);
        acceleration_half(J, 2) = y_dir_sun * G * msun / ((dist_to_sun).^2);

        % iterate through each other planet
        for K = 1:8
            if J ~= K
                % calculate distance to the other planet
                x_dist = position_half(K, 1) - position(J, 1);
                y_dist = position(K, 2) - position(J, 2);
                dist_tot = sqrt((x_dist ^ 2) + (y_dist ^ 2));
                % determine the sign/dir of the other planet
                x_dir = x_dist / dist_tot;
                y_dir = y_dist / dist_tot;
                % acceleration is the accumulated acceleration in that dir plus the new
                acceleration_half(J, 1) = acceleration(J, 1) + x_dir * (G * mass(K)/(dist_tot).^2);
                acceleration_half(J, 2) = acceleration(J, 2) + y_dir * (G * mass(K)/(dist_tot).^2);
            end 
        end 
    end
	
	% position and velocity full step
    for J = 1:8
        % position full step
        position(J,1,I) = position(J,1,I-1) + velocity_half(J,1)*dt;
        position(J,2,I) = position(J,2,I-1) + velocity_half(J,2)*dt;

        % velocity full step
        velocity(J,1) = velocity(J,1,I-1) + acceleration_half(J,1)*dt;
        velocity(J,2) = velocity(J,2,I-1) + acceleration_half(J,2)*dt;
    end 
    
    % Acceleration full step
	% iterate through each planet
	for J = 1:8
		% the new acceleration is initialized as the acceleration toward the sun
		x_dist_sun = xsun-position(J,1,I);
		y_dist_sun = ysun-position(J,2,I);
		dist_to_sun = sqrt(x_dist_sun ^ 2 + y_dist_sun ^ 2);
		x_dir_sun = x_dist_sun / dist_to_sun;
		y_dir_sun = y_dist_sun / dist_to_sun;
		acceleration (J, 1, I) = x_dir_sun * G * msun / ((dist_to_sun).^2);
		acceleration (J, 2, I) = y_dir_sun * G * msun / ((dist_to_sun).^2);
		
        
            % iterate through each other planet
            for K = 1:8
                if J ~= K
                % calculate distance to the other planet
                x_dist = position(K, 1, I) - position(J, 1, I);
                y_dist = position(K, 2, I) - position(J, 2, I);
                dist_tot = sqrt((x_dist ^ 2) + (y_dist ^ 2));
                % determine the sign/dir of the other planet
                x_dir = x_dist / dist_tot;
                y_dir = y_dist / dist_tot;
                % acceleration is the accumulated acceleration in that dir plus the new
                acceleration(J, 1, I) = acceleration(J, 1, I) + x_dir * (G * mass(K)/(dist_tot).^2);
                acceleration(J, 2, I) = acceleration(J, 2, I) + y_dir * (G * mass(K)/(dist_tot).^2);
                end
            end
	end 
end

%% visualization
skip = 1000;

figure(1)
clf
for I = 1:N_steps
    if(mod(I, skip) ~= 0)
        continue
    end
    plot(0,0, 'y. ', 'markersize',30);    
    hold on;
    plot( position(1, 1,I), position(1, 2, I), 'b. ', 'markersize',10);    
    axis([-x08, x08, -x08, x08]); 
    axis equal
    
    
    for P = 2:8    
        plot( position(P, 1,I), position(P, 2, I), 'b. ', 'markersize',10);       
               
    end    
    pause(.01); 
    
end

grid on;
xlabel('t','fontsize',20);
ylabel('x', 'fontsize',20);
zlabel('y', 'fontsize',20);
hold off; 