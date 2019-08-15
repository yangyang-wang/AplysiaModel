% This file computes and plots iPRC. 

xinit=[0.900321164137428;   
    0.083551935956201;
    0.000031666995903;
    0.747647099749367;
    0.246345045901938;
    0.649984712236374;
    -8.273162075117845;
    0.01];

vinit=[0 0 0 0 0 0];

% model=lyttle_model; %creat model
% model.xinit=xinit'; % specify initial value
% T0 = model.findPeriod(4.8,4.9) %find period
T0=4.886087799072266; % Period of LC solution of the Aplysia model when initial value is xinit. 

model = lyttle_model(xinit,vinit,T0);
model.solve; 

% first use find_prc_monodromy to find the initial cond for prc
M = find_prc_monodromy(model);
[V, D]=eig(M');
ind_e1=(abs(diag(D))==max(abs(diag(D)))); %find index where the eigenvalue is 1
zinit=V(:,ind_e1);
dummy=0;
f10=model.lyttle_ODE(dummy, model.yext(end,:),model.checkdomain(model.yext(end,:)));
f10=f10(1:6);
rescale=f10'*zinit;
prcinit=zinit'/(rescale);

model.find_prc(prcinit); % Use find_prc function to compute iPRC
model.plot_prc;          % plot the results

% Below is to estimate relative changes in time in response to static perturbation
% applied to the seaweed: T1=int_{closing time} Z'(t)\part F/\part eps (gamma(t),0)dt = int_{closing time}z_6(t)*(1)

a1=model.yext(:,2);
a2=model.yext(:,3);
ind_close=(a1+a2>0.5);
int=-trapz(model.prct(flipud(ind_close)),model.prc(flipud(ind_close),6));% time is backward,so the integral has the opposite sign 
T1=-int; 
nu=T1/T0;
disp('T1 is')
disp(T1)
disp('nu is')
disp(nu)

%% Plot the results
time_close=model.t(ind_close);
t_clo=time_close(end);

figure
subplot(3,1,1)
plot(model.t,model.yext(:,1:3),'linewidth',2)
hold on
plot([t_clo t_clo], [0 1.1],'g:','linewidth',3)
ylim([-0.01 1.1])

set(gca,'FontSize',18)
legend('a0','a1','a2')
subplot(3,1,2)
plot(model.prct,model.prc(:,1:3),'linewidth',2)
hold on
plot([t_clo t_clo], [-0.01 5.3e+4],'g:','linewidth',3)
ylim([-0.01 5.3e+4])

set(gca,'FontSize',18)
legend('Z_1','Z_2','Z_3')

subplot(3,1,3)
plot(model.prct,model.prc(:,4:6),'linewidth',2)
hold on
plot([t_clo t_clo], [-5 5],'g:','linewidth',3)
ylim([-5 5])
xlim([0 T0])


xlabel('time','interpreter','latex','fontsize',25)
legend('Z_4','Z_5','Z_6')

set(gca,'FontSize',18)