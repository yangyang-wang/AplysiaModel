% Generate Fig.2,Fig.3, Fig.4 (enlargement of Fig.2) and Fig.5

%% Preliminary Simulations
% Find IC at closing and unperturbed period
M=lyttle_model('nu',[]);
M.tmax=15;
M.solve;
xinit = M.y_open_to_close(end,1:8); % start at close
M_period = M.t_open_to_close(end)-M.t_open_to_close(end-1);

% find perturbed solution with applied load fsw+eps
M2=lyttle_model;
eps=0.02;
M2.xinit(8)=M.xinit(8)+eps; M2.tmax=15; M2.solve;

M2.xinit = M2.y_open_to_close(end,1:8);   % start at close
M2.xinit(7) = M.xinit(7); % start with the same seaweed position 
M2_period = M2.t_open_to_close(end)-M2.t_open_to_close(end-1);

M2.tmax = 2*M2_period; % run for two periods
M2.solve

% generate unperturbed solution
vinit=(M2.xinit(1:6)-xinit(1:6))/eps;
M = lyttle_model( 'xinit', xinit, 'vinit', vinit, ...
    'tmax', 2*M_period, 'nu', []);  
M.solve


%estimate linear shift in time T1 and T1_close: need to make eps small to get accurate estimate
T1_num = (M2.t_open_to_close(1)-M.t_open_to_close(1))/eps;
T1_close_num = (M2.t_close_to_open(1)-M.t_close_to_open(1))/eps;


%%  Fig.2

% Fig.2(A,B): Pull out the difference between x and x-perturbed (x2)
% align to time of closing, no rescaling

dt=0.01; % sampling interval -- use for both, to compare
tplot=0:dt:M.tmax;
tplot2=0:dt:M2.tmax;
tidx=find(diff(M.t)>0);
tidx2=find(diff(M2.t)>0);

phi=@(x)-2.598076211353316*x.*(x-1).*(x+1);

% interpolated values of trajectory for comparison

a0i=interp1(M.t(tidx),M.yext(tidx,1),tplot);
a1i=interp1(M.t(tidx),M.yext(tidx,2),tplot);
a2i=interp1(M.t(tidx),M.yext(tidx,3),tplot);
u0i=interp1(M.t(tidx),M.yext(tidx,4),tplot);
u1i=interp1(M.t(tidx),M.yext(tidx,5),tplot);
xpi=interp1(M.t(tidx),M.yext(tidx,6),tplot);
spi=interp1(M.t(tidx),M.yext(tidx,7),tplot);
Fmi=(u0i*M.k0).*phi((M.c0-xpi)/M.w0)+(u1i*M.k1).*phi((M.c1-xpi)/M.w1); % Fmusc
Fmi_pro=(u0i*M.k0).*phi((M.c0-xpi)/M.w0); 
Fmi_re=(u1i*M.k1).*phi((M.c1-xpi)/M.w1); 

a0i2=interp1(M2.t(tidx2),M2.yext(tidx2,1),tplot2);
a1i2=interp1(M2.t(tidx2),M2.yext(tidx2,2),tplot2);
a2i2=interp1(M2.t(tidx2),M2.yext(tidx2,3),tplot2);
u0i2=interp1(M2.t(tidx2),M2.yext(tidx2,4),tplot2);
u1i2=interp1(M2.t(tidx2),M2.yext(tidx2,5),tplot2);
xpi2=interp1(M2.t(tidx2),M2.yext(tidx2,6),tplot2);
spi2=interp1(M2.t(tidx2),M2.yext(tidx2,7),tplot2);
Fmi2=(u0i2*M.k0).*phi((M.c0-xpi2)/M.w0)+(u1i2*M.k1).*phi((M.c1-xpi2)/M.w1); % Fmusc
Fmi_pro2=(u0i2*M.k0).*phi((M.c0-xpi2)/M.w0);
Fmi_re2=(u1i2*M.k1).*phi((M.c1-xpi2)/M.w1);

nx=min(length(xpi),length(xpi2));

%Plot
figure
set(gcf,'Position',[500 800 1200 400])
% Plot trajectories together, aligned from closing
%  plot brain variables
subplot(1,2,1)
plot(tplot(1:nx),a0i(1:nx),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(tplot(1:nx),a1i(1:nx),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot(1:nx),a2i(1:nx), 'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot([M2.t_close_to_open(1) M2.t_close_to_open(1)],[0 1],'m','linewidth',2)
plot([M2.t_close_to_open(2) M2.t_close_to_open(2)],[0 1],'m','linewidth',2)

plot(tplot(1:nx),a0i2(1:nx),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(tplot(1:nx),a1i2(1:nx),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot(1:nx),a2i2(1:nx),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)
legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([tplot(1) tplot(end)])
ylim([0, 1])
M.draw_wall_closing
xlabel('time','interpreter','latex','fontsize',20)


subplot(1,2,2)
plot(tplot(1:nx),u0i(1:nx),'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(tplot(1:nx),u1i(1:nx),'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot(1:nx),xpi(1:nx),'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(tplot(1:nx),Fmi(1:nx),'Color',[0.5 0.5 0.5],'LineWidth',2)

plot([M2.t_close_to_open(1) M2.t_close_to_open(1)],[-0.4 1],'m','linewidth',2)
plot([M2.t_close_to_open(2) M2.t_close_to_open(2)],[-0.4 1],'m','linewidth',2)

plot(tplot(1:nx),u0i2(1:nx),'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(tplot(1:nx),u1i2(1:nx),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot(1:nx),xpi2(1:nx),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(tplot(1:nx),Fmi2(1:nx),'--','Color',[0.5 0.5 0.5],'LineWidth',2)

xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','$F_{\rm musc}$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([tplot(1) tplot(end)])
ylim([-0.4, 1])
M.draw_wall_closing

shg
drawnow

% Fig.2(C,D)
M.plot_var_horizontal
subplot(1,2,1)
ylim([-72 60])
hold on
plot([M2.t_close_to_open(1) M2.t_close_to_open(1)], [-72 60],'m','linewidth',2)
plot([M2.t_close_to_open(2) M2.t_close_to_open(2)], [-72 60],'m','linewidth',2)

% plot([T0_close, T0_close], [-30 30],'b:','linewidth',2)
legend('$a_0$','$a_1$','$a_2$','interpreter','latex','fontsize',15)
title('forward variational ($\mathbf{u}(t)$)','Interpreter','latex','FontWeight','normal','Fontsize',20)            

subplot(1,2,2)
ylim([-5 10])
hold on
plot([M2.t_close_to_open(1) M2.t_close_to_open(1)], [-5 10],'m','linewidth',2)
plot([M2.t_close_to_open(2) M2.t_close_to_open(2)], [-30 30],'m','linewidth',2)
plot([0 t_isrc(end)],[0 0],'k:','linewidth',2)
title('forward variational ($\mathbf{u}(t)$)','Interpreter','latex','FontWeight','normal','Fontsize',20)  
ylim([-8 10])
% add line with predicted initial rate of growth of x2-x1
set(line([0,model.tmax],[0,model.tmax]/model0.br),'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)
legend('$u_0$','$u_1$','$x_r$','$F_{\rm musc}$','interpreter','latex','fontsize',15)


%% Fig.3
figure
plot(tplot(1:nx),Fmi_pro(1:nx),'r','LineWidth',2)
hold on
plot(tplot(1:nx),Fmi_pro2(1:nx),'r:','LineWidth',2)

plot(tplot(1:nx),Fmi_re(1:nx),'b','LineWidth',2)
hold on
plot(tplot(1:nx),Fmi_re2(1:nx),'b:','LineWidth',2)

plot([M2.t_close_to_open(1) M2.t_close_to_open(1)],[-0.6 1],'m','linewidth',2)
plot([M2.t_close_to_open(2) M2.t_close_to_open(2)],[-0.6 1],'m','linewidth',2)

xlabel('time','interpreter','latex','fontsize',20)
ylabel('Muscle Forces','interpreter','latex','fontsize',25)
set(gca,'FontSize',20)
xlim([tplot(1) tplot(end)])
ylim([-0.6, 0.6])
M.draw_wall_closing
text(3,0.5, 'protraction','Interpreter','latex','FontSize',20,'Color','r')
text(3,-0.1, 'retraction','Interpreter','latex','FontSize',20,'Color','b')

%% Fig.5

% Fig.5(A,B): stretch time to run from closing to opening, with piecewise rescaling

tplotc=linspace(0,M.t_close_to_open(1),501); % closed region for first cycle
tploto=linspace(M.t_close_to_open(1),M.t_open_to_close,501); % open region of first cycle

tplotc_2nd=linspace(M.t_open_to_close,M.t_close_to_open(2),501); % closed region of second cycle
tploto_2nd=linspace(M.t_close_to_open(2),max(M.t),501); % open region (to end of sim)

tplot=[tplotc,tploto(2:end),tplotc_2nd(2:end),tploto_2nd(2:end)];

tplot2c=linspace(0,M2.t_close_to_open(1),501); % closed region
tplot2o=linspace(M2.t_close_to_open(1),M2.t_open_to_close(1),501); % open region (to end of sim)

tplot2c_2nd=linspace(M2.t_open_to_close(1),M2.t_close_to_open(2),501); % closed region of second cycle
tplot2o_2nd=linspace(M2.t_close_to_open(2),max(M2.t),501); % open region (to end of sim)


tplot2=[tplot2c,tplot2o(2:end),tplot2c_2nd(2:end),tplot2o_2nd(2:end)];

tidx=find(diff(M.t)>0);
tidx2=find(diff(M2.t)>0);

% interpolated values of trajectory for comparison
a0i=interp1(M.t(tidx),M.yext(tidx,1),tplot);
a1i=interp1(M.t(tidx),M.yext(tidx,2),tplot);
a2i=interp1(M.t(tidx),M.yext(tidx,3),tplot);
u0i=interp1(M.t(tidx),M.yext(tidx,4),tplot);
u1i=interp1(M.t(tidx),M.yext(tidx,5),tplot);
xpi=interp1(M.t(tidx),M.yext(tidx,6),tplot);
spi=interp1(M.t(tidx),M.yext(tidx,7),tplot);


a0i2=interp1(M2.t(tidx2),M2.yext(tidx2,1),tplot2);
a1i2=interp1(M2.t(tidx2),M2.yext(tidx2,2),tplot2);
a2i2=interp1(M2.t(tidx2),M2.yext(tidx2,3),tplot2);
u0i2=interp1(M2.t(tidx2),M2.yext(tidx2,4),tplot2);
u1i2=interp1(M2.t(tidx2),M2.yext(tidx2,5),tplot2);
xpi2=interp1(M2.t(tidx2),M2.yext(tidx2,6),tplot2);
spi2=interp1(M2.t(tidx2),M2.yext(tidx2,7),tplot2);

%
figure
set(gcf,'Position',[500 800 1200 400])

subplot(1,2,1)
plot(tplot,a0i,'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(tplot,a1i,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot,a2i,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(tplot,a0i2,'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(tplot,a1i2,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot,a2i2,'Color',[0.9290 0.6940 0.1250],'LineStyle','--','LineWidth',2)
legend('$a_0$','$a_1$','$a_2$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([tplot(1) tplot(end)])
M.draw_wall_closing
ylim([0, 1])
xlabel('time','interpreter','latex','fontsize',20)

subplot(1,2,2)
plot(tplot,u0i,'Color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(tplot,u1i,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot,xpi,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)

plot(tplot,u0i2,'--','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(tplot,u1i2,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(tplot,xpi2,'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
xlabel('time','interpreter','latex','fontsize',20)
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')
set(gca,'FontSize',20)
xlim([tplot(1) tplot(end)])
M.draw_wall_closing

shg
drawnow

% Fig.5 (C,D)
% Find the iSRC with global uniform or piecewise uniform rescaling; start from close
pt=true;
ICisClose=true;
isUniform=true;
% iSRC with piecewise rescaling, starting from close
[t_isrc_pw_close,x_isrc_pw_close]=lyttle_iSRC(pt,ICisClose,~isUniform);
subplot(1,2,2)
plot([0 tplot(end)],[0 0],'k:','linewidth',2)
ylim([-1 1.5])
legend('$u_0$','$u_1$','$x_r$','Location','southwest','interpreter','latex')


