% Code to analyse the parameters of the Hankel matrix for the Rossler
% attractor
% Written by: Prerna Patil
% Original code reference: Eurika Kaiser
% Date: 25th May 2023

clc
clear
close all

% Matlab plotting parameters
LW = 'linewidth';

% File paths
addpath('../utils');
addpath('../../Data/Rossler/')
figpath = '../MATLABFIGURES/';
datapath = '../DATA/';

%% Load UPO info from the *txt files depending on the case
Case = 'I';
if strcmp(Case,'I')
    filename = 'a043295b2c4CaseI.txt';
    p.a = 0.43295;
    tol = 0.013;
elseif strcmp(Case,'II')
    filename = 'a043295b2c4CaseII.txt';
    p.a = 0.43295;
    tol = 0.013;
elseif strcmp(Case,'III')
    filename = 'a043295b2c4.txt';
    p.a = 0.43295;
    tol = 0.013;
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

%% Plot UPOs
cmap = jet(2*Ncases); cmap = cmap(1:2:end,:);
figure, hold on, box on
counter = 0;
for i = 1:Ncases %(4*6)
    counter = counter + 1;
    subplot(4,3,counter), hold on, box on
    L = PeriodUPO(counter,2);
%     plot3(xdat0(:,1),xdat0(:,2),xdat0(:,3)-4,'-','Color',0.7*ones(1,3),'LineWidth',2)
    plot3(xdatPred(2:L,1,i),xdatPred(2:L,2,i),xdatPred(2:L,3,i),'-','Color',cmap(i,:),'LineWidth',1) 
    view(3)
    title(num2str(i))
    set(gca,'xtick',[], 'ytick', [])
    if counter == (11)
        break
    end
end
set(gcf,'Position',[100 100 1000 600])
set(gcf,'PaperPositionMode','auto')
% Check if the periodic data has been augmented correctly 
% Case 61 for examples has sharp data jumps in the values. The period
% length has not been computed correctly 

%% Sort the UPOs based on the ratio of time spent in y<0 and y>0 
for counter =1:Ncases
    L = 1:PeriodUPO(counter,2);
    % plot(xdat(1:PeriodUPO(counter,2),2,counter));
    tempVec = xdat(1:PeriodUPO(counter,2),2,counter)<0;
    rho = sum(tempVec)/PeriodUPO(counter,2);
    RatioSymDyn(counter,1) = (rho-1)/(rho+1);
   
end
RatioSymDyn(:,2)= 1:Ncases;
% Sort based on the ratio 
RatioSymDyn =  sortrows(RatioSymDyn,1);
UPOSeq = cell(Ncases,1);
for counter=1:Ncases
    UPOSeq{counter} = LegendNames{RatioSymDyn(counter,2)};
end
UPOSeq = [{'chaotic'};UPOSeq];


%% Embedding based on chaotic trajectory
% tau = 10*dt;
tau = 2.5*20*dt;
skip = tau/dt;
% T_height = [0.5, 5, 50, 250];
T_height = 100;
T_width = 1200;
HankelHeight = floor(T_height/tau);
HankelWidth = T_width/tau;
nUPOs = Ncases;
% Consider only the trajectories for the good trajectories represented by
% the vector good_ones
% Can also consider removing the data set (DataJumps) from the original set
colors = jet(2*Ncases); colors = colors(2:2:end,:);
% Set up the video
SaveFigs   = 0;
MakeVideo  = 0;
saveLegend = 0;
if(MakeVideo)
    Filename = sprintf('%s',Case);
    video = VideoWriter(Filename,'MPEG-4');  % Create a VideoWriter object
    video.FrameRate = 20;  % Optional: Set the frame rate
    video.Quality = 100; 
    open(video);
end
for counter = 1:length(HankelHeight)
    close all
    clear H0 H
    % Update embedding
    % Embedding based on chaotic trajectory
    % the value of skip has been hardcoded to 1 
    H0 = CreateHankelMatrix(xdat0(1:skip:end,:),HankelHeight(counter),HankelWidth,1);
    [U0,S0,V0] = fastersvd(H0); % V = S^(-1)*U'*H; % H = U S V'
    
    % Hankel matrix created from data for each UPO
    % The first dimension is made 5 to reduce the memory cost
    VT = zeros(5,HankelWidth,nUPOs);

    for k = 1:nUPOs
        disp(k);
        clear H
        % the value of skip has been hardcoded to 1 
        H(:,:) = CreateHankelMatrix(xdatPred(1:skip:end,:,RatioSymDyn(k,2)),...
            HankelHeight(counter),HankelWidth,1);
        if(counter>1)
            SignChange = SignsH0EVs-sign(V0(1,1:5));
            SignChange = abs(SignChange)*(-1)+1;
            V0(:,1:5) = V0(:,1:5).*SignChange;
            U0(:,1:5) = U0(:,1:5).*SignChange;
        end
        % Embedding coordinates of UPOs based on chaotic embedding modes
        VTemp = S0^(-1)*U0'*H(:,:);
        VT(:,:,k) = VTemp(1:5,:);
    end
    
    if(counter==1)
        % store the signs from the first eigenvectors
        SignsH0EVs = sign(V0(1,1:5));
    end

    MakeAttractorPlot(V0,VT,colors);
%     ColorCodedEmbeddings(V0,VT, UPOColors);
    TitleString = sprintf('T_{height}=%2.2f',T_height(counter));
%     title(TitleString);
%     Filename = sprintf('../MATLABFIGURES/Rossler/UPOs_Case%s/TDC_UPOs_stack%d',Case,EmbedDims(k));
%     saveas(gca,Filename,'epsc');
%     saveas(gca,Filename,'fig');
    hold off
    if(SaveFigs)
        % Save the figure
        Filename = sprintf('RosslerFigures/AllUPOs/Skip%d/%sHeight%d_Width%d',...
            skip,Case,HankelHeight(counter), HankelWidth);
        saveas(gca,Filename,'epsc');
        saveas(gca,Filename,'fig');
    end
    % chaotic trajectory
    if(saveLegend)
        figure
        % chaotic trajectory
        plot3(NaN,NaN,NaN,'-','color',[0.5,0.5,0.5],'LineWidth',3)
        hold on
        for i=1:nUPOs
            plot3(NaN,NaN,NaN,'-','Color',colors(i,:),'LineWidth',3);
        end
        if(strcmp(Case,'I')||strcmp(Case,'II'))
        legend(UPOSeq,'Location','best','box','off','interpreter','tex','NumColumns',1);
        else
            legend(UPOSeq,'Location','best','box','off','interpreter','tex','NumColumns',2);
        end
        xlabel('x'); ylabel('y');
        set(gca,'FontSize',16,'FontName','Times New Roman');
        grid on;
        set(gca, 'visible', 'off');
        set(gcf,'Position',[100 100 150 300])
        Filename = sprintf('RosslerFigures/CaseIIILegend',...
            Case);
        saveas(gca,Filename,'epsc');
        saveas(gca,Filename,'fig');
       
    end
    if(MakeVideo)
        if(counter>1)
            % Set the view of the new figure to match the previous one
            view(az, el);
        end
        for i = 1:3
            % Rotate the camera by 1 degree
            camorbit(1,0,'data',[0,0,1]);
            fname = sprintf('TempFile');
            print('-djpeg','-r200',fname)     % save image with '-r200' resolution
            I = imread([fname '.jpg']);       % read saved image
            frame = im2frame(I);              % convert image to frame
%             writeVideo(wobj,frame); 
            % Capture the frame
%             frame = getframe(gcf);
            % Write the frame to the video
            writeVideo(video,frame);
        end
        % Old view
        [az,el] = view;
    end
    
    disp(['FINI: ',num2str(100*counter/length(HankelHeight)),'%'])
end
% Close the video file
if(MakeVideo)
    close(video);
end
