% This file computes and plots solutions to variational equations. 

xinit = [0.900321164137428;
    0.083551935956201;
    0.000031666995903;
    0.747647099749367;
    0.246345045901938;
    0.649984712236374;
    -8.273162075117845;
    0.01];              % Initial condition for Aplysia model

% model=lyttle_model; %creat model
% model.xinit=xinit'; % specify initial value
% T0 = model.findPeriod(4.8,4.9) %find period

T0=4.886087799072266; % Period of LC solution of the Aplysia model when initial value is xinit. 

vinit_matrix=eye(6);     % initial conditions for the variational equation
M=[];

figure
for i=1:6
    model=lyttle_model(xinit,vinit_matrix(i,:),T0);
    model.solve;
    newcol=(model.yext(end,9:14))';
    M=[M newcol]; %Monodromy matrix
    
    %% Plot output
    subplot(4,2,1)
    plot(model.t,model.yext(:,1:3),'linewidth',2)   % time traces of the three neural variables
    axis([0 T0 -0.1 1.1])
    legend('a_0','a_1','a_2')
    subplot(4,2,2)
    plot3(model.yext(:,1),model.yext(:,2),model.yext(:,3),'linewidth',2) % Projection of solution onto neural phase space (a0,a1,a2)
    axis([-0.1 1.1 -0.1 1.1 -0.1 1.1])
    xlabel('a_0')
    ylabel('a_1')
    zlabel('a_2','rot',0)
    title('Phase plane')
    subplot(4,2,i+2)
    plot(model.t,model.yext(:,9:13),'linewidth',2) % variational dynamics with initial condition given by ith column of vinit_matrix
    xlim([0 T0])
    title(['Column ', num2str(i), ' of \Phi'])
end
disp('The eigenvalues of the Monodromy matrix are')
disp(eig(M))  % Check eigenvalues of the Monodromy matrix: one of the eigenvalues should be 1!
