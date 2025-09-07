function [R2 slope]=Fun_Fred_Reg_Linear(x,y) 

%% Plot the graph %%
plot(x,y,'ko','MarkerFaceColor','k');

%% Compute regression %%
[p,S]=polyfit(x,y,1); 
slope=p(1);p2=p(2); 
% Generate a new x vector with increased resolution and predict corresponding y values %
if min(x)>0
    x_fit=0.8*min(x):abs(max(x)/100):1.25*max(x);
else
    x_fit=1.25*min(x):abs(max(x)/100):0.8*max(x);
end
y_fit=slope*x_fit+p2;
% Find the R² of the regression model %
yresid=y-(slope*x+p2);
SSresid=sum(yresid.^2);
SStotal=(length(y)-1)*var(y);
R2 = round((1 - SSresid/SStotal),3);

%% Plot the regression line on top of the existing graph %%
hold on; plot(x_fit,y_fit,strcat('--','r'));
% Indicate the R² value on the graph %
xpos=min(x_fit)+0.88*(max(x_fit)-min(x_fit)); % x position at 90% of the regression line
y_xpos=y_fit(find(x_fit==max(x_fit))); % y position corresponding to xpos
ypos=y_xpos+0.4*(max(y_fit)-min(y_fit)); % y position 20% higher than y at max x value
text(xpos,ypos,strcat('R²= ' , num2str(round(R2,2))),'FontWeight','bold','FontSize',12,'color','r');
end