% This script uses the Midpoint Algorithm to track N_steps of 8 planets around the sun. 

% The function Find_Acceleration takes a position array and a parameters array 
% with G as the first value and the rest the masses of the planets. It returns an 
% array containing the acceleration of each of the planets in each direction, by
% calculating the gravitational force of each body on each other body and returning
% the summation.

% The positions, velicities, and masses of each body are intialized at the top.

clear all;

% Steps to iterate
N_steps = 100000;
N_skip  = 1000;
num_planets = 9;

% physical constants;
G = 6.67E-11;

% Time array
T = 3.15569e7; % one year in seconds
t = linspace(0,2*T, N_steps);
dt = mean(diff(t));

% Initialze the Sun
xsun = 0; 
ysun = 0;
msun = 1.989*10^30;
vsun = 0;

% Initialize planets, each with an initial x position, y position, and mass

% Mercury
x01  = 57*10^9;
y01  = 0;
vx01 = 0;
vy01 = 47.9*10^3;
m1   = 328.5*10^21;

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

% Initialize arrays with structure array(Planet, Dir, N)
position     = zeros(num_planets ,2,N_steps);
velocity     = zeros(num_planets ,2,N_steps);
acceleration = zeros(num_planets ,2,N_steps);

% Mass array
mass = [m1, m2, m3, m4, m5, m6, m7, m8, msun];

% initialize positions and velocities
position(:, 1, 1) = [x01; x02; x03; x04; x05; x06; x07; x08; xsun];
position(:, 2, 1) = [y01; y02; y03; y04; y05; y06; y07; y08; ysun];
velocity(:, 1, 1) = [vx01; vx02; vx03; vx04; vx05; vx06; vx07; vx08; vsun];
velocity(:, 2, 1) = [vy01; vy02; vy03; vy04; vy05; vy06; vy07; vy08; vsun];

parameters = cat(2, [G], mass);

% initialize accelerations
acceleration(:, :, 1) = Find_Acceleration(position(:, :, 1), parameters);


% Iterate through N_steps
for I = 2: N_steps	

    % Midpoint Algorthm's half step 
    position_half = position(:,:,I-1) + velocity(:,:,I-1)*dt/2;
    velocity_half = velocity(:,:,I-1) + acceleration(:,:,I-1)*dt/2;
    acceleration_half = Find_Acceleration(position_half, parameters);

	
    % Midpoint Algorithm's full step 
    position(:,:,I) = position(:,:,I-1) + velocity_half*dt;
    velocity(:,:,I) = velocity(:,:,I-1) + acceleration_half*dt;
    acceleration(:,:,I) = Find_Acceleration(position(:,:,I), parameters);
end


%% visualization

xmax = max(position(:,1));
xmin = min(position(:,1));
ymax = max(position(:,2));
ymin = min(position(:,2));

xmax = max(abs([xmin,xmax,ymin,ymax]));
xmin = -xmax;
ymax = xmax;
ymin = xmin;


figure(1)
clf
for I = 1:N_steps
    if(mod(I, N_skip) ~= 0)
        continue
    end
    plot(0,0, 'y. ', 'markersize',30);    
    hold on;
    plot( position(1, 1,I), position(1, 2, I), 'b. ', 'markersize',10);    
    axis([-x04, x04, -x04, x04]);     
    
    
    
    for P = 2:4   
        plot( position(P, 1,I), position(P, 2, I), 'b. ', 'markersize',10);       
               
    end    
    pause(.01); 
    axis equal
    
end

grid on;
xlabel('x','fontsize',20);
ylabel('y', 'fontsize',20);
hold off; 