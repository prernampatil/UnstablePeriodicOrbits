function H = MakeAttractorPlot(V0,VT,colors,V1)

Ncases = size(VT,3);

hf=figure
plot3(V0(:,1),V0(:,2),V0(:,3),'-','color',[0.5,0.5,0.5,0.4],'LineWidth',2.1) % chaotic trajectory
hold on 
if(nargin==4)
    plot3(V1(:,1),V1(:,2),V1(:,3),'-','color',[0.5,0.5,0.5,0.4],'LineWidth',2.1) % chaotic trajectory
end
for k = 1:Ncases
    ph(k) = plot3(VT(1,:,k),VT(2,:,k),VT(3,:,k),'-','color',...
        [colors(k,:),1.0],'LineWidth',2);
    
end
% Plotting variables
grid on
axis tight
xlabel('v_1'), ylabel('v_2'), zlabel('v_3')
xticks([-0.02,0,0.02]);
yticks([-0.02,0,0.02]);
zticks([-0.02,0,0.02]);
% xticks([]);
% yticks([]);
% zticks([]);
axis tight
% box off
box on 
% axis off
light               % create a light
lighting gouraud    % preferred method for lighting curved surfaces
%Open figure and keep handle
hf.Color='w'; %Set background color of figure window
set(gca,'FontSize',28,'FontName','Times New Roman')
view([-21,30]);
set(gcf,'Position',[100 100 800 700])
set(gcf,'PaperPositionMode','auto');
end