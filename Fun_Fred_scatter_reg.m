%% The following function imports x and y data and generates a scatter plot with regression line and R² value %%

function Fun_Fred_scatter_reg(x,y,xlbl,ylbl,font)
% Generate dot plot and correlations
plot(x,y,'ko');hold on;
xlabel(xlbl);ylabel(ylbl);
axis equal
set(gca,'TickDir','out');set(gca,'box','off');
set(gca,'FontSize',font);
Fun_Fred_Reg_Linear(x,y);
end