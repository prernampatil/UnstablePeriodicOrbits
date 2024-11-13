% Code to compute the singular values of the Rossler case 
clc
clear
close all

% Matlab plotting parameters
LW = 'linewidth';

% File paths
addpath('../utils');
addpath('../../Data/Rossler/')
figpath = '../MATLABFIGURES/';
datapath = '../../Data/';

%% Load UPO info from the *txt files depending on the case
Case = 'II';
if strcmp(Case,'I')
    filename = 'a0398b2c4.txt';
    p.a = 0.398;
    tol = 0.012;
elseif strcmp(Case,'II')
    filename = 'a043295b2c4CaseI.txt';
    p.a = 0.43295;
    tol = 0.013;
elseif strcmp(Case,'III')
    filename = 'a0492b2c4.txt';
    p.a = 0.492;
elseif strcmp(Case,'IV')
    filename = 'a0556b2c4.txt';
    p.a = 0.556;
end
% import the data based on the filename
DataUPO = importdata(filename);

% Parameter values
p.b = 2; p.c = 4;

% The x coordinate of the initial conditions for the UPOs is taken to be x_
% (one of the fixed points of the attractor)
xP = (p.c-sqrt(p.c^2-4*p.a*p.b))/2;
explanation = {'period', 'number', 'y', 'z', 'word', 'apparition'};
fid = fopen(filename);

% Scan the file for the data (initial conditions of the attractor)
Ncases = textscan(fid, '%f',1); Ncases = Ncases{1};
Ncol = textscan(fid, '%f',1); Ncol = Ncol{1};
if Ncol==6
    DataUPO = fscanf(fid, '%f %f %f %f %f %f\n',[Ncol Ncases]);
elseif Ncol==5
    DataUPO = fscanf(fid, '%f %f %f %f %f\n',[Ncol Ncases]);
end
fclose(fid);

% Consider if the data for any case is zeros
is_not_present = [];
is_not_present_idx = [];
for i = 1:size(DataUPO,2)
    if all(DataUPO(2:end,i)==0)
        is_not_present_idx = [is_not_present_idx,i];
        is_not_present = [is_not_present,DataUPO(1,i)];
    end
end
if isempty(is_not_present)==0
    DataUPO(:,is_not_present) = [];
end
% DataUPO = DataUPO(:,1:5);
Ncases = size(DataUPO,2);

%% Simulate system
% Simulate the Rossler system for chaotic trajectory and the unstable
% periodic orbit trajectories
nvar = 3;
dt   = 0.005;
% Increasing the t beyond 100 causes the periodic orbits to deviate from
% their periodic trajectories (property of UNSTABLE periodic orbits)
t    = 0:dt:70;
N    = length(t);

% The 'RelTol' has been changed from 10^-16 to 10^-13
ROSSLER = @(t, x, p) [-x(2)-x(3);
    x(1) + p.a*x(2);
    p.b + x(3)*(x(1)-p.c)];
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-16*ones(1,nvar));

% Chaotic trajectory
% Obtain the trajectory of the chaotic case
tspan = dt:dt:10000;
x0 = [-2.8585    0.2234    0.2979];  % initial conditions
[~,xdat0]=ode45(@(t,x) ROSSLER(t,x,p),tspan,x0,ode_options);


%% Simulate UPOs
% The UPOs are simulated from their initial conditions read from the data
% files 
% Also create the names of the UPOs 
xdat = zeros(length(t),3,Ncases);
counter = 0;
ICs = zeros(Ncases,3);
LegendNames = cell(Ncases,1);

for i = 1:Ncases
    counter = counter + 1;
    x0    = [xP,DataUPO(3,i),DataUPO(4,i)]';
    [~,xdat(:,:,counter)] = ode45(@(t,x)ROSSLER(t,x,p),t,x0,ode_options);
    disp(['Number of UPOs simulated: ',num2str(counter), ' of ', num2str(Ncases)])
    LegendNames{i} = int2str(DataUPO(5,i));
    % get a list of ICs
    ICs(i,:) = x0;
end

%% Determine the periods
% To perform analysis of the long time embedding of the UPOs the data needs
% to be augmented. This cannot be done by increasing the time span of the
% Rossler simulation as the UPOs deviate from their periodic trajectories. 

% Determine the time period of the UPOs
fs = 1/dt;
PeriodUPO = zeros(Ncases,2);

for counter = 1:Ncases
    
    Traj = xdat(2:end,:,counter);
    IC = xdat(1,:,counter);
    NormDist = sqrt(sum((Traj-IC).^2,2));
    
    % This tolerance value has been determined by looking at the plot of
    % the NormDist of the 4th UPO
    % the tol is based on the case (determined
    I = find(NormDist<tol);  
    
    PeriodUPO(counter,1) = (min(I)-1)*dt;
    PeriodUPO(counter,2) = min(I)-1;
%     PeriodUPO(counter,3) = Iopt;
end

%% Create data augmentation for the periodic orbits 
% Add the periodic orbits and repeat them after one period
L = 1:700;
tspan = 0:dt:10000; % New time units for the periodic orbits
xdatPred = zeros(length(tspan),nvar,Ncases);
for  counter = 1:Ncases
    L = PeriodUPO(counter,2);
    xorbit = xdat(1:L,:,counter);
    torbit = t(1:L); 
    % Use the 1D interpolation function to augment the data 
    tpred = tspan(1:end);
    ypred = [interp1(torbit,xorbit(:,1),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')',...
             interp1(torbit,xorbit(:,2),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')',...
             interp1(torbit,xorbit(:,3),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')'];
    xdatPred(:,:,counter) = ypred;
    disp(['PROGRESS: ',num2str(100*counter/Ncases),'%'])
end    
nUPOs = Ncases;
% Create a hankel of chaotic trajectory
% Update embedding
%% Plotting variables
colors = jet(nUPOs); %colors = colors(1:2:end,:);
%%
clear Si CSSi S0 CSS0
tau = 20*dt;
skip = tau/dt;
T_height = 500;
T_width = 1200;
HankelHeight = T_height/tau;
HankelWidth = T_width/tau;

% Embedding based on chaotic trajectory
H0 = CreateHankelMatrix(xdat0,HankelHeight,HankelWidth,skip);
[~,S0,~] = fastersvd(H0); % V = S^(-1)*U'*H; % H = U S V'
CSS0 = cumsum(diag(S0))/sum(diag(S0));

% Create hankel of UPOs
for i=1:nUPOs
    H = CreateHankelMatrix(xdatPred(:,:,i),HankelHeight,HankelWidth,skip);
    [~,S,~] = fastersvd(H);
    Si(:,i) = diag(S);
    CSSi(:,i) = cumsum(diag(S))/sum(diag(S));

end
% Save the data
% Filename = sprintf("SingularValuesFigs/Rossler%sHeight%d_Width%d.mat",...
%     Case,HankelHeight, HankelWidth);
% save(Filename, "CSS0","S0","CSSi","Si","nUPOs");
% Saved data plot
figure
semilogy(CSS0,'.-','color',[0.5,0.5,0.5,0.4],'MarkerSize',10);% chaotic trajectory);
hold on
for i=1:nUPOs
    semilogy(CSSi(:,i),'.-','MarkerSize',10,'color',[colors(i,:),1.0]);
    hold on
end
set(gca,'FontSize',20,'FontName','Times New Roman');
ylabel('Singular Values of Hankel matrix');
xlim([0 100]);
ylim([10^-2 10^0.1]);
xlabel('Singular Value Index');
ylabel('Normalized cumulative singular value')
% Filename = sprintf("SingularValuesFigs/%sHeight%d_Width%d",...
%     Case, HankelHeight, HankelWidth);
% saveas(gca,Filename,'epsc');
% saveas(gca,Filename,'fig');