% Fig12: Actual and Estimated displacements between perturbed and unperturbed solutions under
% piecewise uniform rescaling. 

%% Find period, open/close time and locations for UNPERTURBED system
model0 = lyttle_model; model0.tmax=20;model0.solve; 
% restart at y_open_to_close
xinit = model0.y_open_to_close(end,1:8);
fsw = model0.xinit(8);
model = lyttle_model('xinit',xinit,'tmax',10); % start at y_open_to_close
model.solve

T0=model.t_open_to_close(1);

y_open_to_close = model.xinit; 
y_close_to_open = model.y_close_to_open(1,:); % position where grasper switches to open    

T0_close = model.t_close_to_open(1); % Compute the total time spent during closing 
T0_open = T0 - T0_close;             % total time spent during open

%% Find nu1's for iSRC

% T1s estimated 
T1 = 8.1225; T1_close=5.1848; 

nu1_close=T1_close/T0_close;      % nu in interior
nu1_open=(T1-T1_close)/T0_open; % compute the nu in boundary

%%  Find IC for iSRC 
eps = 0.001; % fsw -> fsw + eps

% we are free to choose the difference at the open_to_close point to be the IC for var equation
% with perturbation on fsw
model0.xinit(8)=fsw+eps;model0.tmax=10;
model0.solve;

xinit_pert = model0.y_open_to_close(end,1:8);   % this is the initial cond for perturbed LC

% Initial value of iSRC is the difference between open_to_close points
vinit = (xinit_pert(1:6)-xinit(1:6))/eps;  

%% With IC and nu1's, we can compute iSRC.

% Compute the unperturbed solution and the iSRC with piecewise nu's found from local timing response curve
model = lyttle_model( 'xinit', xinit, 'vinit', vinit, ...
    'tmax', T0, 'nu', [nu1_close,nu1_open]);  
model.solve

%% Find period, open/close time and locations for PERTURBED system
% Take open_to_close point as the initial point, compute the perturbed LC over one period
model_pert1 = lyttle_model('xinit',xinit_pert,'tmax',10); % rerun the system to find period
model_pert1.solve;
Teps = model_pert1.t_open_to_close(1);                 % perturbed period

model_pert = lyttle_model('xinit',xinit_pert,'tmax',Teps);
model_pert.solve;

%% Next we separate solutions into two parts, and do interpolation.

% Separate the perturbed solution into two segments, one in close and the other in open

T0_close_pert = model_pert.t_close_to_open(1); % Compute the total time spent during closing 
T0_open_pert = Teps - T0_close_pert;             % total time spent during open

ind_open_pert = (model_pert.t > T0_close_pert); % index for open
time_close_pert=model_pert.t(~ind_open_pert);  % time in close
x_close_pert=model_pert.yext(~ind_open_pert,1:6); % location in close

time_open_pert = [time_close_pert(end); model_pert.t(ind_open_pert)];
x_open_pert=[x_close_pert(end,:); model_pert.yext(ind_open_pert,1:6)];

% Separate the unperturbed solution into two segments, one in close and the other in open
ind_open = (model.t > T0_close);
time_open = model.t(ind_open);
x_open=model.yext(ind_open,1:6);

time_close=[model.t(~ind_open); time_open(1)];  % time in open
x_close=[model.yext(~ind_open,1:6); x_open(1,:)]; % location in open

% Rescale the perturbed time in close to be the same as the unperturbed time for LC in close
[tspan_close, Ind_close] = unique((time_close_pert./T0_close_pert).*T0_close,'stable'); % rescale time such that it has time [0 T0_close]
x_unique_close = x_close_pert(Ind_close,:); % obtain xr value after time is rescaled
x_pert_interp_close = interp1(tspan_close,x_unique_close,time_close); % interpolate using unperturbed solution

% Rescale the perturbed time in open to be the same as the unperturbed
% time for LC in open
[tspan_open, Ind_open] = unique(time_close(end)+ ((time_open_pert-time_close_pert(end))./(time_open_pert(end)-time_close_pert(end))).*T0_open,'stable'); % rescale time such that it has time [T0_close T0]
x_unique_open = x_open_pert(Ind_open,:); % obtain xr value after time is rescaled
x_pert_interp_open = interp1(tspan_open,x_unique_open,time_open); % interpolate using unperturbed solution

% After rescaling, the time LC under perturbation spends in close should be the same as T0_close
T0_close_pert_rescale=tspan_close(end)-tspan_close(1); % T0_close_pert after rescaling

% actual displacement obtained from subtracting unperturbed LC from perturbed LC
displacement=[x_pert_interp_close; x_pert_interp_open ]- [x_close; x_open];

%% plot actual displacement obtained from numerical computation vs approximated disp using iSRC
% model.plot_var

figure
set(gcf,'Position',[0 0 720 520])
subplot(3,1,1)
plot([time_close; time_open],displacement(:,1),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,9),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-1e-3 1.5e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
% ylabel('$a_{0_\varepsilon}(\tau(t))-a_0(t)$','interpreter','latex','fontsize',25)
ylabel('$\Delta a_0(t)$','interpreter','latex','fontsize',25)

legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
title('Piecewise uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(3,1,2)
plot([time_close; time_open],displacement(:,2),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,10),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-6e-3 2e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$\Delta a_1(t)$','interpreter','latex','fontsize',25)
% legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
% title('Piecewise uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(3,1,3)
plot([time_close; time_open],displacement(:,3),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,11),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-5e-4 6e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$\Delta a_2(t)$','interpreter','latex','fontsize',25)
% legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing

figure
set(gcf,'Position',[0 0 720 520])
subplot(3,1,1)
plot([time_close; time_open],displacement(:,4),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,12),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-1e-3 1e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$\Delta u_0(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
title('Piecewise uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(3,1,2)
plot([time_close; time_open],displacement(:,5),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,13),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-5e-4 1e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$\Delta u_1(t)$','interpreter','latex','fontsize',25)
% legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
% title('Piecewise uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(3,1,3)
plot([time_close; time_open],displacement(:,6),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,14),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-1e-3 1.5e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$\Delta x_r(t)$','interpreter','latex','fontsize',25)
% legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing

