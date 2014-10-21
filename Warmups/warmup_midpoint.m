% This script implements the midpoint algorithm to plot the motion of two
% planets. 

% physical constants;
G = 6.67 * (10 ^ -11);

% Planet  1
m1 = 10e24;   % kg
x01 = 0;
y01 = 0;
vx01 = 0;
vy01 = 0;

% Planet 2
m2 = 10e24 ; 
x02 = 10 * 6.38e6; 
y02 = 10 * 6.38e6;
vx02 = 0;  
vy02 = 0;

%% Time Arrays
N_steps = 1000;
t = zeros(1,N_steps);
t(1) = 0; 
dt = 0.01;

%% Planet 1 Arrays

x1 = zeros(1,N_steps);
y1 = zeros(1,N_steps);


vx1 = zeros(1,N_steps);
vy1 = zeros(1,N_steps);

ax1 = zeros(1,N_steps); 
ay1 = zeros(1,N_steps);

x1(1) = x01;
y1(1) = y01;

vx1(1) = vx01;
vy1(1) = vy01;

ax1(1) = ((x02 - x01)/ abs(x02 - x01)) * (G * m2 / sqrt((x01-x02)^2+(y01-y02)^2));
ay1(1) = ((y02 - y01)/ abs(y02 - y01)) * (G * m2 / sqrt((x01-x02)^2+(y01-y02)^2));


%% Planet 2 Arrays

x2 = zeros(1,N_steps);
y2 = zeros(1,N_steps);

vx2 = zeros(1,N_steps);
vy2 = zeros(1,N_steps);

ax2 = zeros(1,N_steps);
ay2 = zeros(1,N_steps);

x2(1) = x02;
y2(1) = y02;

vx2(1) = vx02;
vy2(1) = vy02;

ax2(1) = ((x01 - x02)/ abs(x01 - x02)) * (G * m1 / sqrt((x01-x02)^2+(y01-y02)^2));
ax2(1) = ((y01 - y02)/ abs(y01 - y02)) * (G * m1 / sqrt((x01-x02)^2+(y01-y02)^2));

%% Iterate
for I=2:N_steps

% update time array
t(I) = t(I-1) + dt;

% make a 1/2 step
t_half = t(I-1) + dt/2;

% Planet 1 Half Step
vx_half1 = vx1(I-1) + ax1(I-1)*dt/2;
vy_half1 = vy1(I-1) + ay1(I-1) * dt / 2;
x_half1 = x1(I-1) + vx1(I-1)*dt/2;
y_half1 = y1(I-1) + vy1(I-1)*dt/2;

% Planet 2 Half Step
vx_half2 = vx2(I-1) + ax2(I-1)*dt/2;
vy_half2 = vy2(I-1) + ay2(I-1)*dt/2;
x_half2 = x2(I-1) + vx2(I-1)*dt/2;
y_half2 = y2(I-1) + vy2(I-1)*dt/2;

% Acceleration Half Steps
dist_half = sqrt((x_half1-x_half2)^2+(y_half1-y_half2)^2); 

ax_half1 = ((x_half2 - x_half1)/ abs(x_half2 - x_half1)) * (G * m2 / dist_half);
ay_half1 = ((y_half2 - y_half1)/ abs(y_half2 - y_half1)) * (G * m2 / dist_half);

ax_half2 = ((x_half1 - x_half2)/ abs(x_half1 - x_half2)) * (G * m1 / dist_half);
ay_half2 = ((y_half1 - y_half2)/ abs(y_half1 - y_half2)) * (G * m2 / dist_half);


% Planet 1 Full Step
vx1(I) = vx1(I-1) + ax_half1*dt;
vy1(I) = vy1(I-1) + ay_half1*dt;

x1(I) = x1(I-1) + vx_half1*dt;
y1(I) = y1(I-1) + vy_half1*dt;
    
% Planet 2 Full Step
vx2(I) = vx2(I-1) + ax_half2*dt;
vy2(I) = vy2(I-1) + ay_half2*dt;

x2(I) = x2(I-1) + vx_half2*dt;
y2(I) = y2(I-1) + vy_half2*dt;
  

% Acceleration Full Steps
dist = sqrt((x1(I)-x2(I))^2+(y1(I)-y2(I))^2);
ax1(I) = ((x2(I) - x1(I))/ abs(x2(I) - x1(I))) * (G * m2 / dist);
ay1(I) = ((y2(I) - y1(I))/ abs(y2(I) - y1(I))) * (G * m2 / dist);
ax2(I) =  ((x1(I) - x2(I))/ abs(x1(I) - x2(I))) * (G * m1 / dist);
ay2(I) = ((y1(I) - y2(I))/ abs(y1(I) - y2(I))) * (G * m2 / dist);

    
end
% visualize
figure(1)
clf
hold on;
plot3(t, x1, y1, 'linewidth',2);
plot3(t, x2, y2, 'r', 'linewidth', 2);
axis([0 10 0 10e7 0 10e7])
grid on;
xlabel('t','fontsize',20);
ylabel('x', 'fontsize',20);
zlabel('y', 'fontsize',20);
hold off;
