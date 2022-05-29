% Fig11: Actual and Estimated displacements between perturbed and unperturbed solutions under
% uniform rescaling. 


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
T1 = 8.1225;  % estimated T1
nu1 = T1/T0;

%%  Find IC for iSRC 
eps = 0.02; % fsw -> fsw + eps
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
    'tmax', T0, 'nu', [nu1 nu1]);  
model.solve

%% Find period, open/close time and locations for PERTURBED system
% Take open_to_close point as the initial point, compute the perturbed LC over one period
model_pert1 = lyttle_model('xinit',xinit_pert,'tmax',10); % rerun the system to find period
model_pert1.solve;
Teps = model_pert1.t_open_to_close(1);                 % perturbed period

model_pert = lyttle_model('xinit',xinit_pert,'tmax',Teps);
model_pert.solve;

%% Next we do interpolation.

time = model.t; % time
x=model.yext(:,1:6); % location 

time_pert=model_pert.t;  % perturbed time
x_pert=model_pert.yext(:,1:6); % perturbed location 


% Rescale the perturbed time to be the same as the unperturbed time for LC 
[tspan, Ind] = unique((time_pert./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
x_unique = x_pert(Ind,:); % obtain xr value after time is rescaled
x_pert_interp = interp1(tspan,x_unique,time); % interpolate using unperturbed solution

% After rescaling, the time LC under perturbation spends in close is
ind_close_pert=(x_pert_interp(:,2) + x_pert_interp(:,3) >0.5) & (model.t<4); % index for close
time_close_pert=model.t(ind_close_pert);  % time in close after rescaling
T0_close_pert_rescale=time_close_pert(end)-time_close_pert(1); % T0_close_pert after rescaling

% actual displacement obtained from subtracting unperturbed LC from perturbed LC
displacement=x_pert_interp- x;

%% plot actual displacement obtained from numerical computation vs approximated disp using iSRC
model.plot_var
title('ISRC with uniform rescaling ($\gamma_1$)','Interpreter','latex','FontWeight','normal','Fontsize',20)

figure
subplot(3,1,1)
plot(time,displacement(:,1),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,9),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
% ylim([-1e-3 2e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$a_{0_\varepsilon}(\tau(t))-a_0(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
title('Uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(3,1,2)
plot(time,displacement(:,3),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,11),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
% ylim([-1e-3 2e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$a_{2_\varepsilon}(\tau(t))-a_2(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing

subplot(3,1,3)
plot(time,displacement(:,6),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,14),'r:','linewidth',2)
plot([T0_close_pert_rescale, T0_close_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_close, T0_close], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
% ylim([-1e-3 2e-3])
xlabel('\rm time','interpreter','latex','fontsize',25)
ylabel('$x_{r_\varepsilon}(\tau(t))-x_r(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_closing
