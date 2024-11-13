% Creating the first return map to compute the symbolic dynamics of the
% Rossler attractor 
% Code written by: Prerna Patil 
% Date: 12th February 2024

clc 
clear 
close all 

% Matlab plotting parameters
LW = 'linewidth';

%% Load UPO info from the *txt files depending on the case
p.a = 0.43295;
% Parameter values
p.b = 2; p.c = 4;

% The x coordinate of the initial conditions for the UPOs is taken to be x_
% (one of the fixed points of the attractor)
xP = (p.c-sqrt(p.c^2-4*p.a*p.b))/2;
yP= -(p.c-sqrt(p.c^2-4*p.a*p.b))/(2*p.a);
zP = (p.c-sqrt(p.c^2-4*p.a*p.b))/(2*p.a);
%% Simulate system
% Simulate the Rossler system for chaotic trajectory and the unstable
% periodic orbit trajectories
nvar = 3;
dt   = 0.005;

% The 'RelTol' has been changed from 10^-16 to 10^-13
ROSSLER = @(t, x, p) [-x(2)-x(3);
    x(1) + p.a*x(2);
    p.b + x(3)*(x(1)-p.c)];
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-16*ones(1,nvar));

% Chaotic trajectory
% Obtain the trajectory of the chaotic case
tspan = dt:dt:1000;
x0 = [-2.8585    0.2234    0.2979];  % initial conditions
% x0 = [xP, -3.06, 0.467];
[~,xdat0]=ode45(@(t,x) ROSSLER(t,x,p),tspan,x0,ode_options);

xDiff =[(xdat0(:,1) - xP), xdat0(:,2)]; 
Index=find(xDiff(1:end-1,1)<0 & xDiff(2:end,1) > 0);
%% plot the chaotic trajectory 
figure
plot3(xdat0(:,1),xdat0(:,2),xdat0(:,3),LW,1.5)
xlabel('x'); ylabel('y'); zlabel('z');
hold on 
x = [xP; xP; xP; xP];
y = [yP; -7; -7; xP];
z = [0; 0; 5; 5];
p=fill3(x,y,z,[0.5,0.5,0.5]);
p(1).FaceAlpha = 0.6;
set(gca,'FontName','Times New Roman','FontSize',18);
plot3(xdat0(Index(1:end-1),1),xdat0(Index(1:end-1),2),xdat0(Index(1:end-1),3),'.k','MarkerSize',10);
view([-27, 20]);
grid on 
% filename = sprintf('RosslerPoincareSection');
% saveas(gcf,filename,'epsc');
% saveas(gcf,filename,'fig');
%% Find the intersection of the plane and the Rossler time series 
figure

plot(xdat0(Index(1:end-1),2),xdat0(Index(2:end),2),'.k')
hold on 
plot([-3.09, -3.09],[-6 0],LW,2)
set(gca,'FontName','Times New Roman','FontSize',18);
grid on 
xlabel('y(i)');
ylabel('y(i+1)');
% filename = sprintf('RosslerFRM_a43295');
% saveas(gcf,filename,'epsc');
% saveas(gcf,filename,'fig');

%% Plot the section which divides the symbolic dynamics 
figure
plot3(xdat0(:,1),xdat0(:,2),xdat0(:,3),LW,1.5)
xlabel('x'); ylabel('y'); zlabel('z');
hold on 
x = [-4; 6; 6; -4];
y = [-3.04; -3.04; -3.04; -3.07];
z = [0; 0; 5; 5];
p=fill3(x,y,z,[0.5,0.5,0.5]);
p(1).FaceAlpha = 0.6;
set(gca,'FontName','Times New Roman','FontSize',18);
view([-27, 20]);
grid on 
% filename = sprintf('RosslerSymbolicDivision');
% saveas(gcf,filename,'epsc');
% saveas(gcf,filename,'fig');