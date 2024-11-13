function ColorCodedEmbeddings(V0,VT, UPOColors)

HankelHeight = size(V0,2);
HankelWidth = size(V0,1); 
Ncases = size(VT,3);

FigurePlot = figure('visible','on');

%plot3(V0(:,1),V0(:,2),V0(:,3),'-','color',[0.5,0.5,0.5,0.4],'LineWidth',1.2) % chaotic trajectory
%hold on 
for k = 1:Ncases
    ph(k) = plot3(VT(1,:,k),VT(2,:,k),VT(3,:,k),'-','Color',...
        UPOColors(k,:),'LineWidth',2);
    hold on 
end

% Plotting variables
grid on
axis tight
xlabel('v_1'), ylabel('v_2'), zlabel('v_3')
set(gca,'FontSize',16,'FontName','Times New Roman')
% view([145,15]);
set(gcf,'Position',[100 100 800 700])
set(gcf,'PaperPositionMode','auto');
% Save the figure
% Filename = sprintf('ViswanathFigures/AllUPOSeqLen6Height%d_Width%d',...
%     HankelHeight, HankelWidth);
% saveas(gca,Filename,'epsc');
% saveas(gca,Filename,'fig');

end