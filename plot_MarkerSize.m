function plot_MarkerSize(mu,mu_vec,Y_toPlot,marker,LineForm)
hold on 
plot(mu,Y_toPlot(mu),LineForm,'LineWidth',1.2)
hold on 
plot(mu_vec,Y_toPlot(mu_vec),marker,'LineWidth',1.2)
end