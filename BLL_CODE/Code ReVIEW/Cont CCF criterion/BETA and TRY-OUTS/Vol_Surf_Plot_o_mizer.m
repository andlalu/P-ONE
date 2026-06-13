%% Imp Vol Surface plot-o-mizer

figure
hold all
for t=[1,2,4,5]
surf(std_T_graph,std_K_graph,vols(:,2:end,t)-vols(:,2:end,3),'EdgeColor',((25-t)/25).*[0.85 0.85 0.85],'FaceColor',((25-t)/25).*[0.9 0.9 0.9]);
set(gcf,'Colormap',brewermap(10,'PuBu'))
freezeColors
% surf(std_T_graph,std_K_graph,vols2(:,2:end,t),'EdgeColor',((25-t)/25).*[0.85 0.85 0.85],'FaceColor',((25-t)/25).*[0.9 0.9 0.9]);
% lambda(t)
end
% surf(std_T_graph,std_K_graph,vols(:,2:end,5),'EdgeColor',[.1 .1 .1],'FaceColor',[0 122 192]./255);
% % set(gcf,'Colormap',brewermap(10,'PuBu'))
% freezeColors
% 
% surf(std_T_graph,std_K_graph,vols(:,2:end,18),'EdgeColor',[.8 .8 .8],'FaceColor',[0.95 0.95 .95]);
% set(gcf,'Colormap',brewermap(10,'YlOrRd'))
% colorbar('eastoutside')
% caxis([0 0.5])
% hcb=colorbar;
% set(hcb,'YTick',0.05:0.05:0.5)
% set(hcb,'YTickLabel',{'5%','10%','15%','20%','25%','30%','35%','40%','45%','50%'})

grid on
alpha(.8)
view(-70,15)
set(gca,'yTick',0.8:0.1:1.2)
set(gca,'xTick',[0.1 0.3 0.5 0.7 1])
set(gca,'zTick',-0.1:0.05:0.1)
axis([0.1 1 0.8 1.2 -0.1 0.1])
% % % % % % % % % % % % % set(gca,'yTickLabel',{{'0.8';'OTM'}, '0.9', {'1';'ATM'}, '1.1',{'1.2';'ITM'}})
set(gca,'zTickLabel',{'-10%', '-5%', '0','5%','10%'})
set(findall(gcf,'type','axes'),'fontname','Helvetica','fontsize',12) %'fontweight','bold',
set(findall(gcf,'type','text'),'fontname','Helvetica','fontsize',12)%'fontweight','bold',
title('Impact of \mu_j^Q on option prices')
xlabel('Maturity','rot',40)
xlabh = get(gca,'XLabel');
set(xlabh,'Position',get(xlabh,'Position') + [-0.09 .035 0])
ylabel('Money-ness','rot',-3)
ylabh = get(gca,'YLabel');
set(ylabh,'Position',get(ylabh,'Position') + [0.15 0 0])
zlabel('Implied Volatility');
% h=legend('t=0    , \lambda_t=16.95 (1 jump event)','t=1M, \lambda_t=4.20','t=2M, \lambda_t=1.23','t=3M, \lambda_t=0.53','t=4M, \lambda_t=0.38');
% set(h,'Location','northwest','FontSize',16);
hold off

% % subplot(1,2,2)
% hold all
% surf(std_T_graph,std_K_graph,vols(:,2:end)-vols2(:,2:end),'EdgeColor',[.1 .1 .1]);
% set(gcf,'Colormap',lbmap(2,'Blue'))
% caxis([-0.1 0.1])
% grid on
% colorbar('eastoutside')
% hcb=colorbar;
% alpha(.8)
% set(hcb,'YTick',-0.1:0.1:0.1)
% set(hcb,'YTickLabel',{'-10%','0%','10%'})
% view(-70,15)
% set(gca,'yTick',0.8:0.1:1.2)
% set(gca,'xTick',[0.1 0.3 0.5 0.7 1])
% set(gca,'zTick',-0.1:0.05:0.1)
% axis([0.1 1 0.8 1.2 -0.1 0.1])
% set(gca,'zTickLabel',{'-10%', '-5%', '0%', '5%','10%'})
% set(findall(gcf,'type','axes'),'fontname','Palatino Linotype','fontsize',20) %'fontweight','bold',
% set(findall(gcf,'type','text'),'fontname','Palatino Linotype','fontSize',20)%'fontweight','bold',
% title({'Surface Plot of the Difference in Implied Volatility';'(IV with jump risk - IV with no jump risk)'})
% xlabel('Maturity','rot',40)
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [-0.09 .035 0])
% ylabel('Money-ness','rot',-3)
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [0.15 0 0])
% zlabel({'Implied Volatility';'Difference'});




