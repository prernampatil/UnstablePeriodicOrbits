% Ratio of left to right lobe in time steps
% Code written by: Prerna Patil
% Date: 29th March 2024

clc
clear
close all

LW = 'linewidth';
addpath('../utils/');

%% Simulate Lorenz system (throws away transient)
nVar  = 3;
dt    = 0.001;
t     = dt:dt:27;
% Parameters of the Lorenz attractor
SIGMA = 10;
RHO   = 28;
BETA  = 8/3;
tspan = 0:dt:1000;

LORENZ = @(t, x) [SIGMA*(-x(1)+x(2));
    RHO*x(1)-x(2)-x(1)*x(3);
    -BETA*x(3)+x(1)*x(2)];
ode_options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,nVar));

%% Select the UPO case
% Case = 'Case1';
% Case = 'Case2';
Case = 'Case3';

if(strcmp(Case,'Case1'))
    nUPOs = 23*2;
    foldername = sprintf('../../Data/Lorenz/DATA4/');
    [xdat,PeriodUPO] = GetUPOData(foldername,nUPOs,Case,t,dt,ode_options,LORENZ);
    xdat = xdat(:,:,2:end);
    PeriodUPO = PeriodUPO(2:end,:);
    nUPOs = nUPOs-1;
elseif(strcmp(Case,'Case2'))
    nUPOs= 19;
    foldername = sprintf('../../Data/Lorenz/DATA3/');
    [xdat,PeriodUPO] = GetUPOData(foldername,nUPOs,Case,t,dt,ode_options,LORENZ);

elseif(strcmp(Case,'Case3'))
    nUPOs = 39;
    foldername = sprintf('../../Data/Lorenz/DATA1/');
    [xdat,PeriodUPO] = GetUPOData(foldername,nUPOs,Case,t,dt,ode_options,LORENZ);
end
%% Find the names of the UPOs
UPOSequences = cell(nUPOs,1);
for  counter = 1:nUPOs
    L = PeriodUPO(counter,2);
    xorbit = xdat(1:L,:,counter);
    s = symdyn(xorbit');
    UPOSequences(counter) = {s};
    disp(['PROGRESS: ',num2str(100*counter/nUPOs),'%'])
end

%% Compute the required data
% From the computed data set, it is observed that in case of the Lorenz
% attractor the ratio of the A and B symbols in the UPO is similar to the
% ratio of time spent per period in each lobe by the UPO.
for counter=1:nUPOs
    % find the values of the xdat which are x>0
    tempVec = xdat(1:PeriodUPO(counter,2),1,counter)>0;
    UPOData(counter,1) = sum(tempVec)/(PeriodUPO(counter,2));

    NumAs = length(find(UPOSequences{counter}=='A'));
    NumBs = length(find(UPOSequences{counter}=='B'));
    UPOData(counter,2) = NumBs/(NumBs+NumAs);
    UPOData(counter,3) = abs(UPOData(counter,2) - UPOData(counter,1))/UPOData(counter,2)*100;
end

%% Find the time spent by trajectory in a single loop of one lobe
% This separation is done based on separating the right and left loop
% (crossing the x=0 plane).
% Find the average time per loop for the lorenz attractor
for counter=1:nUPOs
    plot(xdat(1:PeriodUPO(counter,2),1,counter));
    grid on
    hold on 
    tempVec = xdat(1:PeriodUPO(counter,2),1,counter)>0;
    % Find the flip and make a vertical line
    tempVecDiff = abs(diff(tempVec));
    tempIndex = find(tempVecDiff==1);
    for indexCounter = 1: length(tempIndex)
        plot([tempIndex(indexCounter),tempIndex(indexCounter)],[-20,20],'-k',LW,1.5);
    end
    % find the places where derivative is zero
    tempVec = xdat(1:PeriodUPO(counter,2),1,counter);
    tempVecDiff = diff(tempVec)>0;
    tempIndex = find(abs(diff(tempVecDiff))==1);
    for indexCounter = 1: length(tempIndex)
        plot([tempIndex(indexCounter),tempIndex(indexCounter)],[-20,20],'-k',LW,1.5);
    end
    clf
end
close all 
%% Code to evaluate the average values of alpha and beta
AlphaData = zeros(nUPOs,3);
for counter=1:nUPOs
    X = xdat(:,1,counter);
    Xplus = X(X>0);
    Xminus = X(X<0);
    AlphaData(counter,1) = sum(Xplus)/length(Xplus);
    AlphaData(counter,2) = sum(Xminus)/length(Xminus); 
    
end


