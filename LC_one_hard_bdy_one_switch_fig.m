tmax=200;
alpha0=.2;alpha1=.24;
beta0=2;beta1=3;
yinit0=[alpha0,-1];
yinit1=[alpha1,-1];
[t0,y0]=LC_one_hard_bdy_one_switch(yinit0,tmax,alpha0,beta0);
[t1,y1]=LC_one_hard_bdy_one_switch(yinit1,tmax,alpha1,beta1);
figure
set(line([0 0],[-1 2.2]),'color','k','LineStyle','--','LineWidth',2)
hold on
set(line([-1 1.5],[-1 -1]),'color','k','LineWidth',2)
plot(y0(:,1),y0(:,2),'b-','LineWidth',2)
hold on
plot(y1(:,1),y1(:,2),'r-.','LineWidth',2)
plot(alpha0,-1,'b+','MarkerSize',20)
plot(alpha1,-1,'r+','MarkerSize',20)
plot(alpha0,-1,'b.','MarkerSize',20)
plot(alpha1,-1,'r.','MarkerSize',20)
grid on
shg
set(gca,'FontSize',20)
ylim([-1.1 2.2])
legend('Switch','Hard','\gamma_0','\gamma_1')

print -dpdf LC_one_hard_bdy_one_switch_fig.pdf