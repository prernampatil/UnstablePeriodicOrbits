% Code to compare the proximity in embedded and physical space 
% Date: 8th November 2023 

% First run the UnfoldingFoViswanathUPOs.m code to obtain the data of the
% UPOs in the phase space and embedded space. 

% Figure 1 
% Comparison between UPO1 and UPO45
i=1; plot3(xdat(:,1,i),xdat(:,2,i),xdat(:,3,i),'-','color',...
[colors(1,:),1.0],'LineWidth',2);
hold on 
i=33; plot3(xdat(:,1,i),xdat(:,2,i),xdat(:,3,i),'-','color',...
[colors(36,:),1.0],'LineWidth',2); 
view([90 0]) 
xticks([]); yticks([]); zticks([]);
legend(LegendNames{2},LegendNames{34},'Location','southeast' )
set(gca,'FontSize',28,'FontName','Times New Roman');
box on 
% Figure2 
% Comparison between UPO45 and UPO33 
figure;
i=1; plot3(xdat(:,1,i),xdat(:,2,i),xdat(:,3,i),'-','color',...
[colors(1,:),1.0],'LineWidth',2);
hold on 
i=45; plot3(xdat(:,1,i),xdat(:,2,i),xdat(:,3,i),'-','color',...
[colors(45,:),1.0],'LineWidth',2); 
hold on 
i=33; plot3(xdat(:,1,i),xdat(:,2,i),xdat(:,3,i),'-','color',...
[colors(36,:),1.0],'LineWidth',2);
view([90 0]) 
xticks([]); yticks([]); zticks([]);
legend(LegendNames{2},LegendNames{46},LegendNames{34},'Location','southeast' )
set(gca,'FontSize',28,'FontName','Times New Roman');
box on 

% Figure showing proximity in the embedded space 
figure
i=1; plot3(VT(1,:,i),VT(2,:,i),VT(3,:,i),'-','color',...
[colors(1,:),1.0],'LineWidth',2);
hold on 
i=45; plot3(VT(1,:,i),VT(2,:,i),VT(3,:,i),'-','color',...
[colors(45,:),1.0],'LineWidth',2);
i=33; plot3(VT(1,:,i),VT(2,:,i),VT(3,:,i),'-','color',...
[colors(36,:),1.0],'LineWidth',2);
xticks([]); yticks([]); zticks([]);
legend(LegendNames{2},LegendNames{46},LegendNames{34},'Location','southeast' )
set(gca,'FontSize',28,'FontName','Times New Roman');
box on 
xlabel('$V_1$','interpreter','latex'); 
ylabel('$V_2$','interpreter','latex')