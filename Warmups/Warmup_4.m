%% Warm-up 4
% Final project
% position(planet, direction, value)
N_steps = 1000;
position = zeros(2,2,N_steps);
velocity = zeros(2,2,N_steps);
acceleration = zeros(2,2,N_steps);

%% Time Arrays
t = zeros(1,N_steps);
t(1) = 0; 
dt = 0.01;

% physical constants;
G = 6.67 * (10 ^ -11);

%% Planet one
m1 = 10e24; 
x01 = 0;
y01 = 0;

% planet two
m2 = 10e24;
x02 = 10 * 6.38e6;
y02 = 10 * 6.38e6;

% arrays
position(1,1,1) = x01;
position(1,2,1) = y01;

velocity(1,1,1) = 0;
velocity(1,2,1) = 0;

acceleration(1,1,1) = ((x02 - x01)/ abs(x02 - x01)) * (G * m2 / sqrt((x01-x02)^2+(y01-y02)^2));
acceleration(1,2,1) = ((y02 - y01)/ abs(y02 - y01)) * (G * m2 / sqrt((x01-x02)^2+(y01-y02)^2));

% Planet 2 arrays
position(2,1,1) = x02; 
position(2,2,1) = y02; 

velocity(2,1,1) = 0;
velocity(2,2,1) = 0;

acceleration(2,1,1) = ((x01 - x02)/ abs(x01 - x02)) * (G * m1 / sqrt((x01-x02)^2+(y01-y02)^2));
acceleration(2,2,1) = ((y01 - y02)/ abs(y01 - y02)) * (G * m1 / sqrt((x01-x02)^2+(y01-y02)^2));



%% Iterate
for I=2:N_steps

% update time array
t(I) = t(I-1) + dt;

% make a 1/2 step
%t_half = t(I-1) + dt/2;

% Planet 1 Half Step
velocity(1,1,I) = velocity(1,1,I-1) + acceleration(1,1,I-1)*dt;
velocity(1,2,I) = velocity(1,2,I-1) + acceleration(1,2,I-1)*dt;
position(1,1,I) = position(1,1,I-1) + velocity(1,1,I-1)*dt;
position(1,2,I) = position(1,2,I-1) + velocity(1,2,I-1)*dt;
%vx_half1 = vx1(I-1) + ax1(I-1)*dt/2;
%vy_half1 = vy1(I-1) + ay1(I-1) * dt / 2;
%x_half1 = x1(I-1) + vx1(I-1)*dt/2;
%y_half1 = y1(I-1) + vy1(I-1)*dt/2;

% Planet 2 Half Step
velocity(2,1,I) = velocity(2,1,I-1) + acceleration(2,1,I-1)*dt;
velocity(2,2,I) = velocity(2,2,I-1) + acceleration(2,2,I-1)*dt;
position(2,1,I) = position(2,1,I-1) + velocity(2,1,I-1)*dt;
position(2,2,I) = position(2,2,I-1) + velocity(2,2,I-1)*dt;
%vx_half2 = vx2(I-1) + ax2(I-1)*dt/2;
%vy_half2 = vy2(I-1) + ay2(I-1)*dt/2;
%x_half2 = x2(I-1) + vx2(I-1)*dt/2;
%y_half2 = y2(I-1) + vy2(I-1)*dt/2;

% Acceleration Half Steps
%dist_half = sqrt((x_half1-x_half2)^2+(y_half1-y_half2)^2); 

%ax_half1 = ((x_half2 - x_half1)/ abs(x_half2 - x_half1)) * (G * m2 / dist_half);
%ay_half1 = ((y_half2 - y_half1)/ abs(y_half2 - y_half1)) * (G * m2 / dist_half);

%ax_half2 = ((x_half1 - x_half2)/ abs(x_half1 - x_half2)) * (G * m1 / dist_half);
%ay_half2 = ((y_half1 - y_half2)/ abs(y_half1 - y_half2)) * (G * m2 / dist_half);


% Planet 1 Full Step
% vx1(I) = vx1(I-1) + ax_half1*dt;
% vy1(I) = vy1(I-1) + ay_half1*dt;
% 
% x1(I) = x1(I-1) + vx_half1*dt;
% y1(I) = y1(I-1) + vy_half1*dt;
    
% Planet 2 Full Step
% vx2(I) = vx2(I-1) + ax_half2*dt;
% vy2(I) = vy2(I-1) + ay_half2*dt;
% 
% x2(I) = x2(I-1) + vx_half2*dt;
% y2(I) = y2(I-1) + vy_half2*dt;
  

% Acceleration Full Steps
dist = sqrt((position(1,1,I)-position(2,1,I))^2+(position(1,2,I)-position(2,2,I))^2);
acceleration(1,1,I) = ((position(2,1,I)-position(1,1,I))/abs(position(2,1,I)-position(1,1,I)))*(G*m2/dist);
acceleration(1,2,I) = ((position(2,2,I)-position(1,2,I))/abs(position(2,2,I)-position(1,2,I)))*(G*m2/dist);
acceleration(2,1,I) = ((position(1,1,I)-position(2,1,I))/abs(position(1,1,I)-position(2,1,I)))*(G*m1/dist);
acceleration(2,2,I) = ((position(1,2,I)-position(2,2,I))/abs(position(1,2,I)-position(2,2,I)))*(G*m1/dist);
% dist = sqrt((x1(I)-x2(I))^2+(y1(I)-y2(I))^2);
% ax1(I) = ((x2(I) - x1(I))/ abs(x2(I) - x1(I))) * (G * m2 / dist);
% ay1(I) = ((y2(I) - y1(I))/ abs(y2(I) - y1(I))) * (G * m2 / dist);
% ax2(I) =  ((x1(I) - x2(I))/ abs(x1(I) - x2(I))) * (G * m1 / dist);
% ay2(I) = ((y1(I) - y2(I))/ abs(y1(I) - y2(I))) * (G * m2 / dist);

    
end
% visualize
figure(1)
clf
hold on;
x1 = position(1,1,:);
x1 = squeeze(x1);
y1 = position(1,2,:);
y1 = squeeze(y1);
x2 = position(2,1,:);
x2 = squeeze(x2);
y2 = position(2,2,:);
y2 = squeeze(y2);
plot3(t, x1, y1, 'linewidth',2);
plot3(t, x2, y2, 'r', 'linewidth', 2);
axis([0 10 0 10e7 0 10e7])
grid on;
xlabel('t','fontsize',20);
ylabel('x', 'fontsize',20);
zlabel('y', 'fontsize',20);
hold off;
x = position(1,1,:);