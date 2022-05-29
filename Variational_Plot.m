% Compute all variational dynamcis and find the Monodromy matrix
vinit_matrix=eye(6);     % initial conditions for the variational equation
M=[];

model=lyttle_model;
model.solve;

T0 = model.tmax;

figure
set(gcf,'Position',[50 800 800 800])
% Plot output
subplot(4,2,1)
plot(model.t,model.yext(:,1:3),'linewidth',2)   % time traces of the three neural variables
axis([0 T0 -0.1 1.1])
legend('a_0','a_1','a_2')
model.draw_wall_closing

subplot(4,2,2)
plot3(model.yext(:,1),model.yext(:,2),model.yext(:,3),'linewidth',2) % Projection of solution onto neural phase space (a0,a1,a2)
axis([-0.1 1.1 -0.1 1.1 -0.1 1.1])
xlabel('a_0')
ylabel('a_1')
zlabel('a_2','rot',0)
title('Phase plane')
    
for i=1:6
    model=lyttle_model('vinit',vinit_matrix(i,:));
    model.solve;
    newcol=(model.yext(end,9:14))';
    M=[M newcol]; %Monodromy matrix
    
    subplot(4,2,i+2)
    plot(model.t,model.yext(:,9:13),'linewidth',2) % variational dynamics with initial condition given by ith column of vinit_matrix
    xlim([0 T0])
    title(['Column ', num2str(i), ' of \Phi'])
end
disp('The eigenvalues of the Monodromy matrix are')
disp(eig(M))  % Check eigenvalues of the Monodromy matrix: one of the eigenvalues should be 1!


