% Code to visualize the unfolding of the unfolding of the UPOs obtained
% from Viswanath
% Data source: https://dept.math.lsa.umich.edu/~divakar/lorenz/index.html
% Date: 7th November 2023
% Case 1: UPOs of type AnB and ABn
%         DATA4 contains data for 24 periodic orbits
% Case 2: UPOs of type AnBn (Symmetric UPOs)
%         DATA3 contains data for 25 periodic orbits AnBn
% Case 3: UPOs for sequence length less than 7
%         DATA1 contains data about 1375 periodic orbits corresponding to symbol
%         sequences of length 13 of less.

clc
clear
close all
LW = 'linewidth';
addpath('../utils/');

%% Simulate Lorenz system (throws away transient)
nVar  = 3;
dt    = 0.01;
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
Case = 'Case1';
% Case = 'Case2';
% Case = 'Case3';

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

%% Get chaotic trajectory
x0 = [11,14,27];
[~,xdat0]=ode45(@(t,x) LORENZ(t,x),tspan,x0,ode_options);

%% Plotting variables
colors = jet(2*nUPOs); colors = colors(1:2:end,:);
%axis_limits = [-8e-3 8e-3 -10e-3 10e-3 -16e-3 16e-3];

%% Create data augmentation for the periodic orbits
% Add the periodic orbits and repeat them after one period
% tspan = 0:dt:1000; % New time units for the periodic orbits
tspan = 0:dt:1000;
UPOSequences = cell(nUPOs,1);
xdatPred = zeros(length(tspan),nVar,nUPOs);
for  counter = 1:nUPOs
    L = PeriodUPO(counter,2);
    xorbit = xdat(1:L,:,counter);
    s = symdyn(xorbit');
    UPOSequences(counter) = {s};
    xdatPred(:,:,counter) = AugmentSignal(xorbit,L,length(tspan));
    disp(['PROGRESS: ',num2str(100*counter/nUPOs),'%'])
end

%% Create a consise form of the UPO names 
if(strcmp(Case,'Case1')||strcmp(Case,'Case2'))
    LegendNames = cell(nUPOs,1);
    for counter = 1:nUPOs
        NumAs = length(find(UPOSequences{counter}=='A'));
        NumBs = length(find(UPOSequences{counter}=='B'));
        strA = sprintf('A^{%d}',NumAs);
        strB = sprintf('B^{%d}',NumBs);
        LegendNames(counter) = {append(strA,strB)};
        rho = NumAs/NumBs; 
        RatioSymDyn(counter,1) = rho;
    end
    RatioSymDyn(:,2)= 1:nUPOs;
    LegendNames = [{'chaotic'};LegendNames];
    UPOSequences = [{'chaotic'};UPOSequences];

end

%% ------------------------ Data creation completed ---------------- %%
disp('Data augmentation completed for all UPOs');
SaveAugmentedData = 0;
if(SaveAugmentedData)
    save('Data/AugmentedData.mat','xdatPred','PeriodUPO');
    disp('Augmented data is saved in the folder *Data*');
end

%% Sort the UPOs for Case 3 to make the plots prettier 
% Sort the UPOs based on the ratio of A and B 
%% Sort the UPOs based on the ratio of time spent in y<0 and y>0 
if strcmp(Case, 'Case3')
    for counter = 1:nUPOs
        NumAs = length(find(UPOSequences{counter}=='A'));
        NumBs = length(find(UPOSequences{counter}=='B'));
        rho = NumAs/NumBs; 
        RatioSymDyn(counter,1) = rho;
    end
    RatioSymDyn(:,2)= 1:nUPOs;
    % Sort based on the ratio
    RatioSymDyn =  sortrows(RatioSymDyn,1);
    UPOSeq = cell(nUPOs,1);
    for counter=1:nUPOs
        UPOSeq{counter} = UPOSequences{RatioSymDyn(counter,2)};
    end
    UPOSeq = [{'chaotic'};UPOSeq];
    UPOSequences = UPOSeq;
end
%% Create the video for unfolding of the attractor
tau = dt; %[1,2,5,10,20,50]
skip = tau/dt;
T_height = 25;
T_width = 100;
HankelHeight = floor(T_height/dt);
HankelWidth = T_width/tau;
% Set up the video
SaveFigs = 0;
MakeVideo = 0;
saveLegend= 0;
if(MakeVideo)
    Filename = sprintf('%s',Case);
    video = VideoWriter(Filename,'MPEG-4');  % Create a VideoWriter object
    video.FrameRate = 20;  % Optional: Set the frame rate
    video.Quality = 100; 
    open(video);
end
for counter = 1:length(T_height)
    close all
    clear H0 H
    % Update embedding
    % Embedding based on chaotic trajectory
    H0 = CreateHankelMatrix(xdat0,HankelHeight(counter),HankelWidth,skip);
    [U0,S0,V0] = fastersvd(H0); % V = S^(-1)*U'*H; % H = U S V'

    % Hankel matrix created from data for each UPO
    % The first dimension is made 5 to reduce the memory cost
    VT = zeros(5,HankelWidth,nUPOs);

   for k = 1:nUPOs
        disp(k);
        clear H
        H(:,:) = CreateHankelMatrix(xdatPred(:,:,RatioSymDyn(k,2)),...
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

    MakeAttractorPlot(V0,VT,colors);
%     box on 
    TitleString = sprintf('T_{height}=%2.2f',T_height(counter));
%     title(TitleString);
    if(SaveFigs)
        % Save the figure
        Filename = sprintf('Figures/%sHeight%d_Width%d',...
            Case,HankelHeight(counter), HankelWidth);
        saveas(gca,Filename,'epsc');
        saveas(gca,Filename,'fig');
    end
    % chaotic trajectory
    if(saveLegend)
        figure
        % chaotic trajectory
        plot3(NaN,NaN,NaN,'-','color',[0.5,0.5,0.5])
        hold on
        for i=1:nUPOs
            plot3(NaN,NaN,NaN,'-','Color',colors(i,:),'LineWidth',3);
        end
        if(strcmp(Case,'Case1')||strcmp(Case,'Case2'))
        legend(LegendNames,'Location','best','box','off','interpreter','tex','NumColumns',1);
        else
            legend(UPOSequences,'Location','best','box','off','interpreter','tex','NumColumns',1);
        end
        xlabel('x'); ylabel('y');
        set(gca,'FontSize',16,'FontName','Times New Roman');
        grid on;
        set(gca, 'visible', 'off');
        set(gcf,'Position',[100 100 150 700])
        Filename = sprintf('%sLegend',...
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
    disp(['FINI: ',num2str(100*counter/length(T_height)),'%'])
end
% Close the video file
if(MakeVideo)
    close(video);
end