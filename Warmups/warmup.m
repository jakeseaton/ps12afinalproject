%% Final Project Warmup
% Have two objects accelerate toward eachother and collide
clear all;
%% Configuration
%Constants
G = 6.67 * (10 ^ -11);
% arrays
N_steps = 1000;
t = zeros(1,N_steps);
t(1) = 0; 
dt = 25; %0.01;

%% Properties

% Planet  1
m1 = 10e24;   % kg
x01 = 0;
vx01 = 0;

% Mass 2
m2 = 10e24 ; 
x02 = 10 * 6.38e6; 
vx02 = 0;  


%% Planet 1 Arrays

x1 = zeros(1,N_steps);


vx1 = zeros(1,N_steps);


ax1 = zeros(1,N_steps); 

x1(1) = x01;

vx1(1) = vx01;

ax1(1) = ((x02 - x01)/ abs(x02-x01)) * (G * m2 / (x02 - x01) ^ 2);



%% Planet 2 Arrays

x2 = zeros(1,N_steps);

vx2 = zeros(1,N_steps);

ax2 = zeros(1,N_steps);

x2(1) = x02;

vx2(1) = vx02;

ax2(1) = ((x01 - x02)/ abs(x01-x02)) * (G * m1 / (x01 - x02) ^ 2);

%% Evolve the system
for I = 2:N_steps;
    
    % update time array
    t(I) = t(I-1) + dt;
   
    % update velocity arrays
    vx1(I) = vx1(I-1) + ax1(I-1)*dt;
    vx2(I) = vx2(I-1) + ax2(I-1)*dt;
    
    % update position arrays
    x1(I) = x1(I-1) + vx1(I-1)*dt;
    x2(I) = x2(I-1) + vx2(I-1)*dt;
    
    % finally, update acceleration using the current t, velocities, 
    % and positions.
    ax1(I) = ((x2(I) - x1(I))/ abs(x1(I)-x2(I))) * (G * m2 / (x2(I) - x1(I)) ^ 2);
    ax2(I) =  ((x1(I) - x2(I))/ abs(x1(I)-x2(I))) * (G * m1 / (x1(I) - x2(I)) ^ 2);
    
    
end

% visualize
figure(1)
clf
hold on;
plot(t, x1, 'linewidth',2);
plot(t, x2, 'r', 'linewidth', 2);
grid on;
xlabel('t','fontsize',20);
ylabel('x', 'fontsize',20);
hold off;