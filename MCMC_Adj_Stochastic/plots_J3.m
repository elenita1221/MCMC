clear all; 
vars={'J3','u_target','Y3','Y_target','nruns','N','Mx','My','sim_data',...
    'k3','sol3','tol',...
    'E_u_target','E_u_target_square','ts','RelError3'};
load('All_J_N5_10runs-7tol_10x10_noise-1.mat',vars{:});
E_Y3=mean(Y3(k3,:,:),3)
E_Y_target=mean(Y_target(:,:),2)'
E_Y3_square=mean(Y3(k3,:,:).*Y3(k3,:,:),3);
Var_Y3=E_Y3_square-E_Y3.^2
Var_Y_target=mean(Y_target(:,:).*Y_target(:,:),2)'-(E_Y_target).^2

close all;
if ispc, ms = '\'; else ms = '/'; end
ad = cd;
cd(['..',ms,'Porous_Adjoint_Stochastic']);

figure;
%max_k=max([k3,k4,k5]); % maximum no. of iterations for plotting on the x-axis
max_J3=max(J3(:));
xlim([1 k3]); %ylim([0 max_J3]);
plot(1:k3,J3(1:k3),'--g','LineWidth',2);
title(['Convergence of the Cost Functional'])
xlabel('Iterations');ylabel('Cost Functional');
legend('J_3');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig1_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xlim([1 k3]);
plot(1:k3,log10(J3(1:k3)),'--g','LineWidth',2);
title(['Convergence of the Log10(Cost functional)'])
xlabel('Iterations');ylabel('Log10(Cost Functional)');
legend('Log10(J_3)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig2_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for run_counter=1:nruns
 for i=1:Mx
  for j=1:My
    D_target(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y_target(:,run_counter));
    D_Y3(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y3(k3,:,run_counter));
    %D_Y4(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y4(k4,:,run_counter));
    %D_Y5(i,j,run_counter)=feval('Dcoeff',x(i),y(j),N,Y5(k5,:,run_counter));
  end
 end
end
E_D_target=mean(D_target(:,:,:),3);
figure;
surf(sim_data.xcs,sim_data.ycs,E_D_target);title('Mean of the target diffusion coeff D over the runs')
for i=1:Mx
  for j=1:My
    D_E_Y3(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,E_Y3);  
    %D_E_Y4(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y4);
    %D_E_Y5(i,j)=feval('Dcoeff',x(i),y(j),N,E_Y5);
  end
end  
figure;
surf(sim_data.xcs,sim_data.ycs,D_E_Y3(:,:));title('The diffusion coeff D evaluted at the mean of estimated Ys')

figure;
for run_counter=1:nruns
 plot(sim_data.xcs,u_target(:,My/2+1,end,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(sim_data.xcs,E_u_target(:,My/2+1,end),'-r','LineWidth',2);
title('Crossections of the target solution and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig3_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for run_counter=1:nruns
 plot(sim_data.xcs,D_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(sim_data.xcs,E_D_target(:,My/2+1),'-r','LineWidth',2);
title('Crossections of the target diffusion coeff D and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig4_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(sim_data.xcs,E_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,D_E_Y3(:,My/2+1),'--g','LineWidth',2);hold off;
%plot(x,D_E_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
%plot(x,D_E_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Mean of target diffusion coeff versus Mean of estimated diffusion coeff')
legend('E(D_{target})','E(D_3)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig5_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_D_target_square=mean(D_target(:,:,:).*D_target(:,:,:),3);
Var_D_target=E_D_target_square-(E_D_target).^2;
E_D_Y3_square=mean(D_Y3(:,:,:).*D_Y3(:,:,:),3);
Var_D_Y3=E_D_Y3_square-(D_E_Y3).^2;
% E_D_Y4_square=mean(D_Y4(:,:,:).*D_Y4(:,:,:),3);
% Var_D_Y4=E_D_Y4_square-(D_E_Y4).^2;
% E_D_Y5_square=mean(D_Y5(:,:,:).*D_Y5(:,:,:),3);
% Var_D_Y5=E_D_Y5_square-(D_E_Y5).^2;
figure;
plot(sim_data.xcs,Var_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,Var_D_Y3(:,My/2+1),'--g','LineWidth',2);hold off;
% plot(x,Var_D_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
% plot(x,Var_D_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Variance of target diffusion coeff versus Variance of estimated diffusion coeff')
legend('Var(D_{target})','Var(D_3)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig6_',num2str(Mx)]],'-depsc')

z=repmat(E_Y3(:),1,1000)+sqrt(Var_Y3(:))*randn(1,1000);
figure;
for i=1:5 subplot(1,5,i);[N,X]=hist(z(i,:));bar(X,N/sum(N));set(gca,'XLim',[0 1]);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 sol3_matrix=zeros(Mx,My,length(ts),nruns);%sol4_matrix=zeros(Mx+2,My+2,nruns);sol5_matrix=zeros(Mx+2,My+2,nruns);
 for run_counter=1:nruns
   sol3_matrix(:,:,:,run_counter)=sol3{run_counter};
%   sol4_matrix(:,:,run_counter)=sol4{run_counter};
%   sol5_matrix(:,:,run_counter)=sol5{run_counter};
 end
 E_sol3_matrix=mean(sol3_matrix(:,:,end,:),4);
% E_sol3_square=mean(sol3_matrix(:,:,:).*sol3_matrix(:,:,:),3);
% E_sol4_matrix=mean(sol4_matrix(:,:,:),3);%already calculated as E_u4
% E_sol5_matrix=mean(sol5_matrix(:,:,:),3);%already calculated as E_u5
 figure;
 surf(sim_data.xcs,sim_data.ycs,E_sol3_matrix(:,:));title('Mean of the estimated solution over the runs')
 figure;
 surf(sim_data.xcs,sim_data.ycs,E_u_target(:,:,end));title('Mean of the target solution over the runs')
 figure;
 plot(sim_data.xcs,E_u_target(:,My/2+1,end),'-r','LineWidth',2);hold on;
 plot(sim_data.xcs,E_sol3_matrix(:,My/2+1,end),'--g','LineWidth',2);hold on;
% plot(x,E_sol4_matrix(:,My/2+1),'-.b','LineWidth',2);hold on;
% plot(x,E_sol5_matrix(:,My/2+1),':m','LineWidth',2);hold off;
% title('Crossections of Mean of target solution versus Mean of estimated solution ')
% legend('E(u_{target})','E(u_3)','E(u_4)','E(u_5)');
% print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig7_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Var_u_target=E_u_target_square-(E_u_target).^2;  %Variance of u_target
% Var_u3=E_sol3_square-(E_sol3_matrix).^2;  %Variance of solution u3
% Var_u4=E_u4_square-(E_u4).^2;  %Variance of solution u4
% Var_u5=E_u5_square-(E_u5).^2;  %Variance of solution u5
% figure;
% plot(x,Var_u_target(:,My/2+1),'-r','LineWidth',2);hold on;
% plot(x,Var_u3(:,My/2+1),'--g','LineWidth',2);hold on;
% plot(x,Var_u4(:,My/2+1),'-.b','LineWidth',2);hold on;
% plot(x,Var_u5(:,My/2+1),':m','LineWidth',2);hold off;
% title('Crossections of Variance of target solution versus Variance of estimated solution ')
% legend('Var(u_{target})','Var(u_3)','Var(u_4)','Var(u_5)');
% print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig8_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E_f_target=mean(f_target(:,:,:),3);
% figure;
% for run_counter=1:nruns
%  plot(x,f_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
%  hold on;
% end
% plot(x,E_f_target(:,My/2+1),'-r','LineWidth',2);
% title('Crossections of the target forcing function versus its mean')
% hold off;
% print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol'],ms,['fig9_',num2str(Mx)]],'-depsc')