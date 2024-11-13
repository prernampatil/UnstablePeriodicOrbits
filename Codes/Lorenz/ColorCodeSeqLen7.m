% Color coding sequences less than sequence lenght 7 
% Sort UPOs based on the time spend in A side or B side
% Written by: Prerna Patil
% Date: 9th Oct 2023

clc
clear
close all
LW = 'linewidth';
addpath('../utils');
%% Simulate Lorenz system (throws away transient)
nVar  = 3;
dt    = 0.01;
t     = dt:dt:25;
% Parameters of the Lorenz attractor
SIGMA = 10;
RHO   = 28;
BETA  = 8/3;
tspan = 0:dt:1000;
% N     = length(tspan);

LORENZ = @(t, x) [SIGMA*(-x(1)+x(2));
    RHO*x(1)-x(2)-x(1)*x(3);
    -BETA*x(3)+x(1)*x(2)];
ode_options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,nVar));

k=0;
nUPO = 21;
% Obtain the length of period for each UPO
PeriodUPO = zeros(nUPO,2);

for counter=1:nUPO
    k=k+1
    filename = sprintf('../../Data/Lorenz/DATA1/orbit%d.dat',counter);
    orbit=load(filename);
    x0 = orbit(:,1);
    [~,xdat(:,:,k)] = ode45(@(t,x)LORENZ(t,x),t,x0,ode_options);
    filename = sprintf('../../Data/Lorenz//DATA1/T%d.dat',counter);
    PeriodUPO(k,1) =  load(filename);
    PeriodUPO(k,2) = ceil(PeriodUPO(k,1)/dt);
end
% Obtain the symbolic sequences for each UPO
UPOSequences = cell(nUPO,1);
for  counter = 1:nUPO
    L = PeriodUPO(counter,2);
    xorbit = xdat(1:L,:,counter);
    s = symdyn(xorbit');
    UPOSequences(counter) = {s};
end

%% Create a chaotic trajectory
x0 = [10,20,25];
[~,xdat0]=ode45(@(t,x) LORENZ(t,x),tspan,x0,ode_options);


%% Plotting variables
colors = jet(2*nUPO); colors = colors(1:2:end,:);
colorblue = colors(end,:);
colorred = colors(1,:);
colorneutral = [0.5,0.5,0.5];

% Determine the colors for each UPO based on the number of the A and Bs in
% the symbolic sequence
UPOColors = zeros(nUPO,3);
for counter=1:nUPO
    tempstr = UPOSequences{counter};
    num_A = sum(tempstr=='A');
    num_B = sum(tempstr=='B');
    if(num_A>num_B)
        UPOColors(counter,:) = colorblue;
    elseif(num_A<num_B)
        UPOColors(counter,:) = colorred;
    elseif(num_A==num_B)
        UPOColors(counter,:) = colorneutral;
    end
end

%% Create data augmentation for the periodic orbits
% Add the periodic orbits and repeat them after one period
% tspan = 0:dt:1000; % New time units for the periodic orbits
tspan = 0:dt:1000;
xdatPred = zeros(length(tspan),nVar,nUPO);
for  counter = 1:nUPO
    L = PeriodUPO(counter,2);
    torbit = t(1:L);
    xorbit = xdat(1:L,:,counter);
    % Use the 1D interpolation function to augment the data
    tpred = tspan(1:end);
    ypred = [interp1(torbit,xorbit(:,1),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')',...
        interp1(torbit,xorbit(:,2),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')',...
        interp1(torbit,xorbit(:,3),mod(tpred,PeriodUPO(counter,1)),'pchip','extrap')'];
    xdatPred(:,:,counter) = ypred;
    disp(['PROGRESS: ',num2str(100*counter/nUPO),'%'])
end

%% Sort the UPOs based on period length
PeriodUPO(:,3) = 1:nUPO;
PeriodUPO = sortrows(PeriodUPO,2);
SortedUPO = PeriodUPO(:,3);
UPOSequences = UPOSequences(SortedUPO);
UPOColors  = UPOColors(SortedUPO,:);

%% Plot the results for different values of T_height
% xdatPred = xdatPred(1:2:end,:,:);

tau = 2*dt;
skip = tau/dt;
T_height = 15;
T_width = 100;
HankelHeight = floor(T_height/dt);
HankelWidth = T_width/tau;
% Set up the video
MakeVideo = 0;
saveLegend= 0;
if(MakeVideo)
    video = VideoWriter('a.avi');  % Create a VideoWriter object
    video.FrameRate = 15;  % Optional: Set the frame rate
    open(video);
end
for counter = 1:length(T_height)
    close all
    clear H0 H
    % Update embedding
    % Embedding based on chaotic trajectory
    %     H0 = getHankelMatrix(xdat0,HankelHeight(counter),skip);
    H0 = CreateHankelMatrix(xdat0,HankelHeight(counter),HankelWidth,skip);
    [U0,S0,V0] = fastersvd(H0); % V = S^(-1)*U'*H; % H = U S V'

    % Hankel matrix created from data for each UPO
    % Compute coordinates
    % H = zeros(HankelHeight(counter),HankelWidth);
    % VT = zeros(HankelHeight(counter),HankelWidth,Ncases);
    % The first dimension is made 5 to reduce the memory cost
    VT = zeros(5,HankelWidth,nUPO);

    for k = 1:nUPO
        k
        clear H
        H(:,:) = CreateHankelMatrix(xdatPred(:,:,SortedUPO(k)),...
            HankelHeight(counter),HankelWidth,skip);
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

    ColorCodedEmbeddings(V0,VT, UPOColors);
    TitleString = sprintf('T_{height}=%2.2f',T_height(counter));
    title(TitleString);
    % Chaotic trajectory
    if(saveLegend)
        h = zeros(3,1);
        h(1)=plot3(NaN,NaN,NaN,'-','Color',colorblue,'LineWidth',3);
        h(2)=plot3(NaN,NaN,NaN,'-','Color',colorred,'LineWidth',3);
        h(3)=plot3(NaN,NaN,NaN,'-','Color',colorneutral,'LineWidth',3);
        legend(h,{'AnB','BnA','AnBn'},'Location','EastOutside','box','off')

        %         Filename = sprintf('SomeFig',...
        %             HankelHeight, HankelWidth);
        %         saveas(gca,Filename,'epsc');
        %         saveas(gca,Filename,'fig');
    end
    if(MakeVideo)
        if(counter>1)
            % Set the view of the new figure to match the previous one
            view(az, el);
        end

        for i = 1:2
            % Rotate the camera by 1 degree
            camorbit(1,0,'data',[0,0,1]);

            % Capture the frame
            frame = getframe(gcf);

            % Write the frame to the video
            writeVideo(video,frame);
        end
        % Old view
        [az,el] = view;
    end
    %     MakeAttractorPlotProjections(V0,VT,colors);

    disp(['FINI: ',num2str(100*counter/length(T_height)),'%'])
end
view(-145, 1.86);