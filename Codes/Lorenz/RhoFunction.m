% Code to show the collection of UPOs in the ratio of symbols
% Written by: Prerna Patil

clc
close all
clear

colors = jet(2*40); colors = colors(1:2:end,:);
colorblue = colors(1,:);
colorred = colors(end,:);
colorneutral = [0.5,0.5,0.5];

LW = 'linewidth';
addpath('../utils/');


%% Figure 1 
rho = 1:0.01:20;
rhofunc = (rho-1)./(rho+1);

plot(rho, rhofunc,'-k',LW,2.5);
hold on
set(gca,'FontSize',20,'FontName','Times New Roman');
xlabel('$\rho$','Interpreter','latex');
ylabel('$(\rho-1)/(\rho+1)$','Interpreter','latex');
grid on 

for i=1:20
    x = [0,i];
    y = [(i-1)/(i+1),(i-1)/(i+1)];
    plot(x,y,'-',LW,1,'color',colorred);

    x = [i,i];
    y = [-1,(i-1)/(i+1)];
    plot(x,y,'-',LW,1,'color',colorred);
end
% Save figure 
% Filename = sprintf('Rho_1to20');
%         saveas(gca,Filename,'epsc');
%         saveas(gca,Filename,'fig');
%% Figure 2 
figure 
rho = 1/200:0.0001:1;

rhofunc = (rho-1)./(rho+1);

plot(rho, rhofunc,'-k',LW,2.5);
hold on
set(gca,'FontSize',20,'FontName','Times New Roman');
xlabel('$\rho$','Interpreter','latex');
ylabel('$(\rho-1)/(\rho+1)$','Interpreter','latex');
ylim([-1 1]);
grid on 
for i=1:20
    x = [0,1/i];
    y = [(1/i-1)/(1/i+1),(1/i-1)/(1/i+1)];
    plot(x,y,'-',LW,1,'color',colorblue);

    x = [1/i,1/i];
    y = [-1,(1/i-1)/(1/i+1)];
    plot(x,y,'-',LW,1,'color',colorblue);

end
% Save figure 
% Filename = sprintf('Rho_frac1to20');
%         saveas(gca,Filename,'epsc');
%         saveas(gca,Filename,'fig');
