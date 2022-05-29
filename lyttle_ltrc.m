% This file computes and plots local timing response curve starting from close
% created July 2, 2020


function T1_int = lyttle_ltrc(eps1,pl)

% pl: plot if ture. 

% Consider two regions, above and below the wedge
% Generate local timing response curve for the region I that is above the wedge, usng the function find_prc in LC_in_square
% and compute the piecewise relative change in frequency nu_above_wedge and nu_below_wedge
if nargin<1
    eps1 = 0.0001;
end
if nargin<2
    pl = true;
end

% For different eps values, find y_open_to_close, y_close_to_open and the total closing time. 
if eps1 == 0.001
% Need to run model0 for the transients to decay for eps1 other than 0.002
% and to find y_open_to_close as the initial condition
model0 = lyttle_model('eps1',eps1,'tmax',20); 
model0.solve;
% restart at y_open_to_close
model = lyttle_model('eps1',eps1,'tmax',5,'xinit',model0.y_open_to_close(end,1:8)); % start at y_open_to_close
model.solve

elseif eps1 == 0.0001
% with default eps1 = 0.0001
model0 = lyttle_model; % start at y_open_to_close
model0.solve;
model = lyttle_model('xinit',model0.y_open_to_close(end,1:8)); % start at y_open_to_close
model.solve;
end

y_open_to_close = model.xinit; 
y_close_to_open = model.y_close_to_open(1,:); % position where grasper switches to open            
T0_close = model.t_close_to_open(1); % Compute the total time spent during closing 

% To only compute the lTRC in the closing, compute the solution from
% x_open_to_close to x_close_to_open
model = lyttle_model('xinit',y_open_to_close,'tmax',T0_close,'eps1',eps1);
model.solve;

% To compute the lTRC over the full cycle, compute the solution with IC at x_liftoff on [0,T0]
% model = stick_slip('xinit', x_liftoff, 'tmax',T0);
% model.solve;

% Compute the boundary value of lTRC at the close_to_open point, lTRCinit,
% This will be the initial condition for the lTRC since we will integrate the adjoint equation backward in time
lTRCinit0=[0,1,1,0,0,0]'; % the normal vector to open/close wall: y+z=0.5
dummy=0;
f10=model.lyttle_ODE(dummy,y_close_to_open',0); % vector field at the close_to_open point
f10 = f10(1:6);
rescale=f10'*lTRCinit0;
lTRCinit=-lTRCinit0'/(rescale); % the value for the lTRC at the close_to_open point 

model.find_prc(lTRCinit); 

%% Compute T1_close
% T1_int = 0;
% % Claim lTRC(xin)* (xin_pert-xin)/eps = 0, this is the first term in T1_int where in
% % % denotes open-to-close since perturbation during open vanishes so
% % trajectory converges back to the same open-to-close point:
% % x^eps_open_to_close = x^0_open_to_close
% % 
% % eps = 0.001; % FSW -> FSW+eps
T1_int_1=0;
% % 
% integral of the lTRC over the closing region,
% there is a minus sign in front of the integral since time is backward
int=-trapz(model.prct,model.prc(:,6)/model.br);% time is backward,so the integral has the opposite sign 

T1_int=int; 
T1_int = T1_int_1 + T1_int; % relative shift in time in region I
% nu1_int=T1_int/T0_int; % relative change in frequency in the interior
% 
% disp('nu1_int is')
% disp(nu1_int)
% 
if pl
%%
% model.plot_prc;          % visualize the lTRC result

figure
subplot(2,1,1)
plot(model.prct,model.prc(:,1:3),'linewidth',2)
hold on
plot([model.t_exit_a0 model.t_exit_a0],[-20000 2],'k:','linewidth',2)
xlabel('$\rm time$','interpreter','latex','fontsize',30)
% ylabel('$\eta^{\rm close}$','interpreter','latex','fontsize',30,'rot',0)
legend('$\eta_{a_0}$','$\eta_{a_1}$','$\eta_{a_2}$','interpreter','latex','fontsize',25)
xlim([0 model.tmax])
ylim([-20000, 0])
model.draw_wall_closing
title('$\eta^{\rm close}$','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)

subplot(2,1,2)
plot(model.prct,model.prc(:,4:6),'linewidth',2)
hold on
plot([model.t_exit_a0 model.t_exit_a0],[-3.5 2.1],'k:','linewidth',2)
xlabel('$\rm time$','interpreter','latex','fontsize',30)
% ylabel('$\eta^{\rm close}$','interpreter','latex','fontsize',30,'rot',0)
legend('$\eta_{u_0}$','$\eta_{u_1}$','$\eta_{x_r}$','interpreter','latex','fontsize',25)
xlim([0 model.tmax])
ylim([-3.5 2.1])
set(gca,'FontSize',18)
model.draw_wall_closing

end
end