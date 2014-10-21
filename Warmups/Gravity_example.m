%% Two Dimensional Example for N-body simulations
% RSK, 4/16/14
% This script numerically integrates the acceleration a 3-D system of 
% N objects. 

% extract basic parameters
num_steps = 5000; % number of time steps
N = 2;            % number of objects
skip = 10;        % number of frames to skip

% intitial conditions
G    = 6.67E-11;  % N m^2/kg^2, Newton's constant
AU   = 149597871*1e3; % m, one astronomical unit, i.e. distance from earth to sun
Msun = 1.9891E30; % kg, mass of the sun
yr   = 3.15569E7; % s, 

% object 1
m1    = 1; % kg
x1_0  = 1*AU;  % km
y1_0  = 0*AU;  % km

x2_0  = -1*AU; % km
y2_0  = 0*AU;  % km

vx1_0  = 0;    % km/s
vy1_0  = 30e3; % km

vx2_0  = 0;     % km
vy2_0  = -30e3; % km

% allocate memory t array using 1-D array
t = linspace(0,.5*yr, num_steps);
dt = mean(diff(t));

% pre-allocate memory for x,v,a using 2-D array
% first index indicates which object
% 2nd index gives the time step
x = zeros(N,num_steps);
y = zeros(N,num_steps);

vx = zeros(N,num_steps);
vy = zeros(N,num_steps);

ax = zeros(N,num_steps);
ay = zeros(N,num_steps);

% mass array
m  = zeros(N,1);

% initialize mass arrays
m(1) = Msun;
m(2) = Msun;

% initalize position arrays
x(:,1) = [-AU; AU];
y(:,1) = [ 0; 0];  

% initalize the velocity arrays
vx(:,1)= [0; 0];
vy(:,1)= [3e4; -3e4];  

% m/s^2  initialize the acceleration arrays
d = sqrt((x(1,1) -x(2,1))^2 + (y(1,1) - y(2,1))^2);
ax(1,1) = G*m(2)*(x(2,1)-x(1,1))/d^3;
ax(2,1) = G*m(1)*(x(1,1)-x(2,1))/d^3;

ay(1,1) = G*m(2)*(y(2,1)-y(1,1))/d^3;
ay(2,1) = G*m(1)*(y(1,1)-y(2,1))/d^3;


for n = 1:num_steps      
    
    % look ahead, use slope at midpoint of interval
    t_half = t(n) + dt/2;
    
    x_half = x(:,n) + vx(:,n)*dt/2;
    y_half = y(:,n) + vy(:,n)*dt/2;
    
    vx_half = vx(:,n) + ax(:,n)*dt/2;
    vy_half = vy(:,n) + ay(:,n)*dt/2;
    
    % compute accel at midpoint
    d_half = sqrt((x_half(1) - x_half(2))^2 + ...
                  (y_half(1) - y_half(2))^2);
    ax_half = [G*m(2)*(x_half(2) - x_half(1))/d^3; 
               G*m(1)*(x_half(1) - x_half(2))/d^3];
    
    ay_half = [G*m(2)*(y_half(2) - y_half(1))/d^3; 
               G*m(1)*(y_half(1) - y_half(2))/d^3];
    
    
    % using the previous midpoint values, compute the full time step        
    x(:,n+1) = x(:,n) + vx_half*dt;
    y(:,n+1) = y(:,n) + vy_half*dt;
    
    vx(:,n+1) = vx(:,n) + ax_half*dt;
    vy(:,n+1) = vy(:,n) + ay_half*dt;
    
    % compute new accel
    d = sqrt((x(1,n+1) - x(2,n+1))^2 +(y(1,n+1) - y(2,n+1))^2);
    ax(:,n+1) = [G*m(2)*(x(2,n+1) - x(1,n+1))/d^3;
    	         G*m(1)*(x(1,n+1) - x(2,n+1))/d^3];

    ay(:,n+1) = [G*m(2)*(y(2,n+1) - y(1,n+1))/d^3;
    	         G*m(1)*(y(1,n+1) - y(2,n+1))/d^3];    
end


%% visualize

figure(1)
xmin = 1.2*min([min(x(:)),min(y(:))]);
xmax = 1.2*max([max(x(:)),max(y(:))]);
ymin = xmin;
ymax = xmax;

for n = 1:length(t)
    
    % skip frames to speed up animation
    if( mod(n,skip) == 0 && n > 1) 
        pause(.01);
        continue
    end
        
    % plot current positions
    plot(x(1,n), y(1,n), 'b. ', 'markersize',30);
    axis([xmin, xmax, ymin, ymax]);
    hold on;
    plot(x(2,n), y(2,n), 'r. ', 'markersize',30);
    
    if( n == 1)
        pause();        
    end
    hold off;
end