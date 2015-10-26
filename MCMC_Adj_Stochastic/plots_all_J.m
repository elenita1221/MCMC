clear all; % Plots 9 figures
vars={'J3','J4_new','J5_new','u_target','Y3','Y4','Y5','Y_target','nruns','N','Mx','My','sim_data',...
    'k3','k4','k5','sol3','sol4','sol5','tol','E_u3','E_u3_square','E_u4','E_u4_square','E_u5','E_u5_square',...
    'E_u_target','E_u_target_square','ts','quant3','quant4','quant5','noise'};
load('All_J_N5_50runs-7tol_10x10_noise-3.mat',vars{:});
E_Y3=mean(Y3(k3,:,:),3)
E_Y4=mean(Y4(k4,:,:),3)
E_Y5=mean(Y5(k5,:,:),3)
E_Y_target=mean(Y_target(:,:),2)'
E_Y3_square=mean(Y3(k3,:,:).*Y3(k3,:,:),3);
E_Y4_square=mean(Y4(k4,:,:).*Y4(k4,:,:),3);
E_Y5_square=mean(Y5(k5,:,:).*Y5(k5,:,:),3);
Var_Y3=std(Y3(k3,:,:),0,3).^2; %Var_Y3=E_Y3_square-E_Y3.^2
Var_Y4=std(Y4(k4,:,:),0,3).^2; %Var_Y4=E_Y4_square-E_Y4.^2
Var_Y5=std(Y5(k5,:,:),0,3).^2; %Var_Y5=E_Y5_square-E_Y5.^2
Var_Y_target=std(Y_target(:,:),0,2).^2; %Var_Y_target=mean(Y_target(:,:).*Y_target(:,:),2)'-(E_Y_target).^2
std_Y3=std(Y3(k3,:,:),0,3)
std_Y4=std(Y4(k4,:,:),0,3)
std_Y5=std(Y5(k5,:,:),0,3)
std_Y_target=std(Y_target(:,:),0,2)'

close all;
if ispc, ms = '\'; else ms = '/'; end
ad = cd;
cd(['..',ms,'MCMC_Adj_Stochastic']);

figure;
axis tight %sets the axis limits to the range of the data
max_k=max([k3,k4,k5]); % maximum no. of iterations for plotting on the x-axis
max_J3=max(J3(:));max_J4=max(J4_new(:));max_J5=max(J5_new(:));max_J=max([max_J3,max_J4,max_J5]);
xlim([1 max_k]);
plot(1:k3,log10(J3(1:k3)),'--g','LineWidth',2);hold on;
plot(1:k4,log10(J4_new(1:k4)),'-.b','LineWidth',2);hold on;
plot(1:k5,log10(J5_new(1:k5)),':m','LineWidth',2);hold off;
title(['Convergence of the Log10(Cost functional)'])
xlabel('Iterations');ylabel('Log10(Cost Functional)');
legend('Log10(J_3)','Log10(J_4)','Log10(J_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig1_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
xlim([1 max_k]); ylim([0 max_J]);
plot(1:k3,J3(1:k3),'--g','LineWidth',2);hold on;
plot(1:k4,J4_new(1:k4),'-.b','LineWidth',2);hold on;
plot(1:k5,J5_new(1:k5),':m','LineWidth',2);hold off;
title(['Convergence of the Cost Functional'])
xlabel('Iterations');ylabel('Cost Functional');
legend('J_3','J_4','J_5');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig2_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
set(gca,'YScale','log');

max_quant3=max(quant3(:));max_quant4=max(quant4(:));max_quant5=max(quant5(:));
max_quant=max([max_quant3,max_quant4,max_quant5]);
min_quant3=min(quant3(:));min_quant4=min(quant4(:));min_quant5=min(quant5(:));
min_quant=min([min_quant3,min_quant4,min_quant5]);
xlim([1 max_k]);% 
plot(1:k3,quant3(1:k3),'--g','LineWidth',2);hold on;
plot(1:k4,quant4(1:k4),'-.b','LineWidth',2);hold on;
plot(1:k5,quant5(1:k5),':m','LineWidth',2);hold off;
title(['Mean Convergence in L^2 norm of the estimated solution'])
xlabel('Iterations');ylabel('Quantity of interest ');
legend('\int||E(u_3)-E(u_{target})||^2dt','\int||E(u_4)-E(u_{target})||^2dt','\int||E(u_5)-E(u_{target})||^2dt');
set(gca,'YTick',[min_quant, max_quant]);
mticks = get(gca,'YTick');
for mtc = 1:length(mticks)
   mticks2{mtc} = sprintf('%7.4e',mticks(mtc));
end
set(gca,'YTickLabel',mticks2);
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig3_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for run_counter=1:nruns
 for i=1:Mx
  for j=1:My
    D_target(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y_target(:,run_counter));
    D_Y3(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y3(k3,:,run_counter));
    D_Y4(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y4(k4,:,run_counter));
    D_Y5(i,j,run_counter)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y5(k5,:,run_counter));
  end
 end
end
E_D_target=mean(D_target(:,:,:),3);
figure;
surf(sim_data.xcs,sim_data.ycs,E_D_target);title('Mean of the target diffusion coeff D over the runs')
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig4_',num2str(Mx)]],'-depsc')
for i=1:Mx
  for j=1:My
    D_E_Y3(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,E_Y3);  
    D_E_Y4(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,E_Y4);
    D_E_Y5(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,E_Y5);
  end
end  
%figure;
%surf(sim_data.xcs,sim_data.ycs,D_E_Y3(:,:));title('The diffusion coeff D evaluted at the mean of estimated Ys')

figure;
for run_counter=1:nruns
 plot(sim_data.xcs,u_target(:,My/2+1,end,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(sim_data.xcs,E_u_target(:,My/2+1,end),'-r','LineWidth',2);
title('Crossections of the target solution and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig5_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
for run_counter=1:nruns
 plot(sim_data.xcs,D_target(:,My/2+1,run_counter),'-.b','LineWidth',2);
 hold on;
end
plot(sim_data.xcs,E_D_target(:,My/2+1),'-r','LineWidth',2);
title('Crossections of the target diffusion coeff D and its mean')
hold off;
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig6_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(sim_data.xcs,E_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,D_E_Y3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(sim_data.xcs,D_E_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(sim_data.xcs,D_E_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Mean of target diffusion coeff versus Mean of estimated diffusion coeff')
legend('E(D_{target})','E(D_3)','E(D_4)','E(D_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig7_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_D_target_square=mean(D_target(:,:,:).*D_target(:,:,:),3);
Var_D_target=E_D_target_square-(E_D_target).^2;
E_D_Y3_square=mean(D_Y3(:,:,:).*D_Y3(:,:,:),3);
Var_D_Y3=E_D_Y3_square-(D_E_Y3).^2;
E_D_Y4_square=mean(D_Y4(:,:,:).*D_Y4(:,:,:),3);
Var_D_Y4=E_D_Y4_square-(D_E_Y4).^2;
E_D_Y5_square=mean(D_Y5(:,:,:).*D_Y5(:,:,:),3);
Var_D_Y5=E_D_Y5_square-(D_E_Y5).^2;
figure;
plot(sim_data.xcs,Var_D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,Var_D_Y3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(sim_data.xcs,Var_D_Y4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(sim_data.xcs,Var_D_Y5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Variance of target diffusion coeff versus Variance of estimated diffusion coeff')
legend('Var(D_{target})','Var(D_3)','Var(D_4)','Var(D_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig8_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol3_matrix=zeros(Mx,My,length(ts),nruns);sol4_matrix=zeros(Mx,My,length(ts),nruns);sol5_matrix=zeros(Mx,My,length(ts),nruns);
for run_counter=1:nruns
  sol3_matrix(:,:,:,run_counter)=sol3{run_counter};
  sol4_matrix(:,:,:,run_counter)=sol4{run_counter};
  sol5_matrix(:,:,:,run_counter)=sol5{run_counter};
end
E_sol3_matrix=mean(sol3_matrix(:,:,end,:),4);
E_sol3_square=mean(sol3_matrix(:,:,end,:).*sol3_matrix(:,:,end,:),4);
%E_sol4_matrix=mean(sol4_matrix(:,:,:,:),4);%already calculated as E_u4
%E_sol5_matrix=mean(sol5_matrix(:,:,:,:),4);%already calculated as E_u5
figure;
plot(sim_data.xcs,E_u_target(:,My/2+1,end),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,E_u3(:,My/2+1,end),'--g','LineWidth',2);hold on;
plot(sim_data.xcs,E_u4(:,My/2+1,end),'-.b','LineWidth',2);hold on;
plot(sim_data.xcs,E_u5(:,My/2+1,end),':m','LineWidth',2);hold off;
title('Crossections of Mean of target solution versus Mean of estimated solution ')
legend('E(u_{target})','E(u_3)','E(u_4)','E(u_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig9_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Var_u_target=E_u_target_square(:,:,end)-(E_u_target(:,:,end)).^2;  %Variance of u_target
Var_u3=E_sol3_square-(E_sol3_matrix).^2;  %Variance of solution u3
Var_u4=E_u4_square(:,:,end)-(E_u4(:,:,end)).^2;  %Variance of solution u4
Var_u5=E_u5_square(:,:,end)-(E_u5(:,:,end)).^2;  %Variance of solution u5
figure;
plot(sim_data.xcs,Var_u_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,Var_u3(:,My/2+1),'--g','LineWidth',2);hold on;
plot(sim_data.xcs,Var_u4(:,My/2+1),'-.b','LineWidth',2);hold on;
plot(sim_data.xcs,Var_u5(:,My/2+1),':m','LineWidth',2);hold off;
title('Crossections of Variance of target solution versus Variance of estimated solution ')
legend('Var(u_{target})','Var(u_3)','Var(u_4)','Var(u_5)');
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig10_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
%subplot(2,2,1);
axes('Position',[0.13 0.5838 0.3347 0.3412]); 
surf(sim_data.xcs,sim_data.ycs,E_u_target(:,:,end));title('Mean of target solution ');zu=max(zlim);zlim([0 zu]);
%get(gca,'Position') % get the position of the current axis
%subplot(2,2,2);
axes('Position',[0.5703 0.5838 0.3347 0.3412]);
surf(sim_data.xcs,sim_data.ycs,E_sol3_matrix(:,:));title('Mean of estimated solution using J_3');zu=max(zlim);zlim([0 zu]);
%get(gca,'Position') % get the position of the current axis
%subplot(2,2,3);
axes('Position',[0.13 0.11 0.3347 0.3412]);
surf(sim_data.xcs,sim_data.ycs,E_u4(:,:,end));title('Mean of estimated solution using J_4');zu=max(zlim);zlim([0 zu]);
%get(gca,'Position') % get the position of the current axis
%subplot(2,2,4);
axes('Position',[0.5703 0.11 0.3347 0.3412]);
surf(sim_data.xcs,sim_data.ycs,E_u5(:,:,end));title('Mean of estimated solution using J_5');zu=max(zlim);zlim([0 zu]);
%get(gca,'Position') % get the position of the current axis
print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig11_',num2str(Mx)]],'-depsc')
clear zlim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z=repmat(E_Y_target(:),1,1000)+sqrt(Var_Y_target(:))*randn(1,1000);
figure;
for i=1:5 
   subplot(1,5,i);
   [N,X]=hist(z(i,:));
   bar(X,N/sum(N));
   set(gca,'XLim',[0.2 0.8]); 
   set(gca,'FontSize',8);
end
title('Histogram for Y_{target}');
% axes('Position',[0.13 0.11 0.1237 0.815]); 
% [N,X]=hist(z(1,:));bar(X,N/sum(N));%set(gca,'XLim',[0 1]);
% axes('Position',[0.2928 0.11 0.1237 0.815]); 
% [N,X]=hist(z(2,:));bar(X,N/sum(N));%set(gca,'XLim',[0 1]);
% axes('Position',[0.4556 0.11 0.1237 0.815]); 
% [N,X]=hist(z(3,:));bar(X,N/sum(N));%set(gca,'XLim',[0 1]);
% axes('Position',[0.6184 0.11 0.1237 0.815]); 
% [N,X]=hist(z(4,:));bar(X,N/sum(N));%set(gca,'XLim',[0 1]);
% axes('Position',[0.7812 0.11 0.1237 0.815]); 
% [N,X]=hist(z(5,:));bar(X,N/sum(N));%set(gca,'XLim',[0 1]);
% print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig12_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z3=repmat(E_Y3(:),1,1000)+sqrt(Var_Y3(:))*randn(1,1000);
figure;
for i=1:5 
   subplot(1,5,i);
   [N,X]=hist(z3(i,:));
   bar(X,N/sum(N));
   set(gca,'XLim',[0.2 0.8]); 
   set(gca,'FontSize',8);
end
title('Histogram for Ys using J_3')
%print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig13_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z4=repmat(E_Y4(:),1,1000)+sqrt(Var_Y4(:))*randn(1,1000);
figure;
for i=1:5 
   subplot(1,5,i);
   [N,X]=hist(z4(i,:));
   bar(X,N/sum(N));
   set(gca,'XLim',[0.2 0.8]);
   set(gca,'FontSize',8);
end
title('Histogram for Ys using J_4')
%print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig14_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z5=repmat(E_Y5(:),1,1000)+sqrt(Var_Y5(:))*randn(1,1000);
figure;
for i=1:5 
   subplot(1,5,i);
   [N,X]=hist(z5(i,:));
   bar(X,N/sum(N));
   set(gca,'XLim',[0.2 0.8]);
   set(gca,'FontSize',8);
end
title('Histogram for Ys using J_5')
%print(gcf,[['Figures_All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol','_noise',num2str(log10(noise))],ms,['fig15_',num2str(Mx)]],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%