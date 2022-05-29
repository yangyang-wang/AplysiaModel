% Fig.8: Use the piecewise iSRC with Piecewise Uniform Rescaling [nu1_close,
% nu1_open] to compute robustness delta_x1/delta_x0 - T1/T0 

% compare the analytic results with numerical results 


%% Find period, close time and seaweed intake for UNPERTURBED system
model1 = lyttle_model; model1.tmax=20; model1.solve; 

% restart at y_open_to_close
xinit = model1.y_open_to_close(end,1:8);
fsw = model1.xinit(8);
model2 = lyttle_model('xinit',xinit,'tmax',10); % start at y_open_to_close
model2.solve

T0=model2.t_open_to_close(1); % find the total period
T0_close = model2.t_close_to_open(1); % find the total time spent during closing 

% find grasper positions at open-to-close and at close-to-open
xr_open_to_close = model2.yext(1,6);
xr_close_to_open = model2.y_close_to_open(1,6);

% intake of seaweed before perturbation
delta_x0 = -(xr_close_to_open - xr_open_to_close); 

%% Find T1, T1_close, and nu1_close from iPRC and lTRC. nu1_close is need for computing iSRC
T1=8.0778; % computed from T1=lyttle_iPRC; this is the linear shift in full period, should approach T1 = (Teps-T0)/eps as eps->0
T1_close=5.1817; % computed from T1_close = lyttle_ltrc; this is the linear shift in time spent in close

nu1_close=T1_close/T0_close;      % nu1 in close

% nu1_close = T1/T0; % if using global nu1 to compute iSRC, then the estimated robustness doesn't agree well with the numerical results. 

%%  Find IC for iSRC (at open-to-close location) = (open-to-close_eps - open-to-close_0)/eps where eps is perturbation on fsw
eps=0.0001; 
model0=lyttle_model;
model0.xinit(8)=fsw+eps; % apply perturbation eps to fsw
model0.tmax=10;
model0.solve;

% Initial value of iSRC is the difference between open_to_close points
vinit = (model0.y_open_to_close(end,1:6)-xinit(1:6))/eps;  

%% With IC and nu1's, we can compute iSRC with piecewise nu's found from local timing response curve

% Since we only need the trajectory over the closed region to compute robustness,
% we just compute traj for [0 T0_close] for which nu is nu1_close

model = lyttle_model( 'xinit', xinit, 'vinit', vinit, ...
    'tmax', T0_close, 'nu', nu1_close);  
model.solve

% The iSRC along xr direction over the closed region is
gamma1_xr_close = model.yext(:,14); 

% Relative shift in delta_x0 can be estimated by isrc as 
delta_x1 = -(gamma1_xr_close(end) - gamma1_xr_close(1)); 

%% compute relative change in S, task fitness, for different eps values

% eps_range = (0.0001:0.0001:0.001); % fsw -> fsw + eps
eps_range = linspace(0.0001,0.01,30);
RCS_num_part1 = zeros(size(eps_range));
RCS_num_part2 = zeros(size(eps_range));
RCS_ana_part1 = zeros(size(eps_range));
RCS_ana_part2 = zeros(size(eps_range));
RCS_num = zeros(size(eps_range));
RCS_ana = zeros(size(eps_range));

model_pert0 = lyttle_model; % construct a model; in the loop we will vary eps

for i = 1:length(eps_range)
    
    eps = eps_range(i);

%% Run perturbed system isomg open_to_close point as the initial point, find Teps, and grasper position at open-to-close and close-to-open
model_pert0.xinit(8)=fsw+eps;
model_pert0.tmax=10;
model_pert0.solve;

xinit_pert = model_pert0.y_open_to_close(end,1:8);   % Take open-to-close as the initial cond for perturbed LC

model_pert = lyttle_model('xinit',xinit_pert,'tmax',10); % rerun the system to find period
model_pert.solve;

Teps = model_pert.t_open_to_close(1);                 % perturbed period

% grasper positions at open-to-close and at close-to-open
xr_open_to_close_pert = model_pert.xinit(1,6);
xr_close_to_open_pert = model_pert.y_close_to_open(1,6);

%% now compute relative change in S task fitness (RCS): (S(eps)-S0)/S0 = eps*delta_x1/delta_x0 - eps*T1/T0

delta_xeps = -(xr_close_to_open_pert - xr_open_to_close_pert); % intake of seaweed after perturbation

RCS_num_part1(i) = (delta_xeps-delta_x0)/delta_x0; % compute part 1 numerically: (delta_xeps-delta_x0)/delta_x0
RCS_ana_part1(i) = eps*delta_x1/delta_x0; % estimated value for robustness_part1: eps*delta_x1/delta_x0

RCS_num_part2(i) = (Teps-T0)/T0; % compute part 2 numerically: (Teps-T0)/T0
RCS_ana_part2(i) = eps*T1/T0; % estimated value for robustness_part1: eps*T1/T0

RCS_num(i) = RCS_num_part1(i) - RCS_num_part2(i); % robustness computed numerically
RCS_ana(i) = RCS_ana_part1(i) - RCS_ana_part2(i); % robustness estimated from iSRC and iPRC

% disp('numerical robustness is')
% disp(robustness_num)
% disp('semi-analytic robustness is')
% disp(robustness_ana)
end

% the robustness ((S(eps)-S0)/S0)/eps is
robustness_ana = RCS_ana./eps_range;

disp('delta_x1/delta_x0 is') % this is delta_x1/delta_x0
disp(delta_x1/delta_x0)

disp('T1/T0 is') % this is delta_x1/delta_x0
disp(T1/T0)

disp('estimated robustness is')
disp(robustness_ana(1))

%% plot actual displacement obtained from numerical computation vs approximated disp using iSRC

% plot robustness (numeric and analytic) w.r.t. eps
figure
plot(eps_range,RCS_num,'b*','MarkerSize',5,'linewidth',2)
hold on
plot(eps_range,RCS_ana,'rO','MarkerSize',10,'linewidth',2)
xlim([eps_range(1) eps_range(end)])
ylabel('relative reduction of $S$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
text(0.001,-0.01, 'fsw $\to$ fsw + $\varepsilon$' ,'FontWeight','normal','interpreter','latex','fontsize',20)
xlabel('$\varepsilon$','interpreter','latex','fontsize',25)
% ylim([-1e-2 0])

% plot robustness (numeric and analytic) w.r.t. eps/fsw
figure
plot(eps_range/fsw,RCS_num,'b*','MarkerSize',5,'linewidth',2)
hold on
plot(eps_range/fsw,RCS_ana,'rO','MarkerSize',10,'linewidth',2)

xlabel('relative increase of mechanical load ($\varepsilon/\rm fsw$)','interpreter','latex','fontsize',25)
ylabel('relative reduction of $S$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
% ylim([-1e-2 0])