% eps3=0.00046; fsw=0.045; S0=0.096983;
% eps3=2e-5; fsw=0.045; S0=0.093518;



figure; 
plot(model.t,model.yext(:,6),'b','linewidth',2);
hold on; 
plot([0,model.tmax],[model.s1 model.s1],'k:','linewidth',2)
plot([0,model.tmax],[model.s3 model.s3],'k:','linewidth',2)
model.draw_wall_closing
xlim([0 6])
set(gca,'FontSize',20)
xlabel('$\rm time$','interpreter','latex','fontsize',30)
ylabel('$x_r$','interpreter','latex','fontsize',30,'rot',0)