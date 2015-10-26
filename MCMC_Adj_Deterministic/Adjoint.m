clear all %Adjoint Parabolic Code for the Deterministic Case
N=5; RelError(1)=1e1; eps=10000; beta=1e-6; tol=1e-16; eps=2*eps/3;
Mx=40; My=40; %%Omega=[-0.05, 0.05, -0.035, 0.035];grid size
%u_ext='u_exact'; f_src='f_source'; scheme='ec_model2'; 
Y_target(1:N) = 0.5; %L=1;L is just a constant in the formula for u_exact
mySeed = 10; rng(mySeed); 
k=1;Y(k,1:N) =0.6+rand(1,N)/100;  %initial guess for 1st iteration
f=zeros(Mx,My); ts = linspace(0,0.02,20); dt=ts(2)-ts(1);
sim_data.nx = Mx; sim_data.ny = My;  
sim_data.xl = -0.05; sim_data.xr = 0.05;
sim_data.dx = (sim_data.xr-sim_data.xl)./sim_data.nx;
sim_data.yl = -0.035; sim_data.yr = 0.035; 
sim_data.dy = (sim_data.yr-sim_data.yl)./sim_data.ny;
sim_data.xes =linspace(sim_data.xl,sim_data.xr,sim_data.nx+1);
sim_data.yes =linspace(sim_data.yl,sim_data.yr,sim_data.ny+1);
sim_data.xcs = (sim_data.xes(1:end-1)+sim_data.xes(2:end))./2;
sim_data.ycs = (sim_data.yes(1:end-1)+sim_data.yes(2:end))./2;
[xM,yM] = meshgrid(sim_data.xcs,sim_data.ycs);
ecini = double((xM.^2+yM.^2)<0.01^2);
[u_target,~,~]= ec_model2(N,Y_target,ecini,ts,sim_data,f,'plot_soln',0);
%[x,y,dx,dy,u,int_k,params]=feval(scheme,N,ts,L,Y(k,1:N),Y_target,Omega,Mx,My,u_ext,f_src);%SE elliptic solver

% for m=1:(Mx+2)
% for n=1:(My+2)
%    u_target(m,n,1:length(ts))=feval(u_ext,x(m),y(n),N,ts,L,Y_target);
% end
% end
[u,int_D,params]= ec_model2(N,Y(k,1:N),ecini,ts,sim_data,f,'plot_soln',0);

for tc=1:length(ts)
  J(k,tc) =sum(sum(((u(:,:,tc)-u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2; %cost functional
end
J_new(k)=dt*sum(J(k,1:length(ts)))+beta*int_D/2;

for i=1:Mx
 for j=1:My
    cos_prod(i,j,1:N)=cos((1:N)*pi*sim_data.xcs(i)).*cos((1:N)*pi*sim_data.ycs(j))./(1e+3);
    cos_sum(i,j)=sum(cos_prod(i,j,:))./(1e+3);
 end
end 

while RelError(k)>tol
  eps=3*eps/2; 
  k=k+1;
  [eta,~,~]= ec_model2_adjoint(N,Y(k-1,1:N),ecini,ts,sim_data,u-u_target,'plot_soln',0);
  %[~,~,~,~,eta,~,~]=feval('parabolic5pt_adjoint',N,ts,L,Y(k-1,1:N),u,u_target,Omega,Mx,My,u-u_target);%AE elliptic solver
  for tc=1:length(ts)
     gradx_u = [diff(u(:,:,tc),1,1)./sim_data.dx];
     grady_u = [diff(u(:,:,tc),1,2)./sim_data.dy];
     gradx_eta=[diff(eta(:,:,tc),1,1)./sim_data.dx];
     grady_eta=[diff(eta(:,:,tc),1,2)./sim_data.dy];
     matrix1=gradx_u(:,:).*gradx_eta(:,:);
     matrix2=grady_u(:,:).*grady_eta(:,:);
     temp1=-sum(sum(cos_sum(1:Mx-1,:).*matrix1(:,:)))*sim_data.dx*sim_data.dy
     temp2=-sum(sum(cos_sum(:,1:My-1).*matrix2(:,:)))*sim_data.dx*sim_data.dy
     dJdY(k-1,tc)=dt*(temp1+temp2);
  end 
  for i=1:N
   dJdY_new(k-1,i)=sum(dJdY(k-1,1:length(ts)))+beta*sum(sum(params.D(:,:).*cos_prod(:,:,i)))*sim_data.dx*sim_data.dy;
  end
  Y(k,1:N)=Y(k-1,1:N)-eps*dJdY_new(k-1,1:N);
  [u,int_D,params]= ec_model2(N,Y(k,1:N),ecini,ts,sim_data,f,'plot_soln',0);
  %[x,y,dx,dy,u,int_D,params]=feval(scheme,N,ts,L,Y(k,1:N),Y_target,Omega,Mx,My,u_ext,f_src);
  for tc=1:length(ts)
    J(k,tc) =sum(sum(((u(:,:,tc)-u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2;
  end
  J_new(k)=dt*sum(J(k,1:length(ts)))+beta*int_D/2;
  
  while J_new(k)>=J_new(k-1)
     eps=eps/10
     if eps<1e-15 disp('Algorithm stagnated');break; end; 
     Y(k,1:N)=Y(k-1,1:N)-eps*dJdY_new(k-1,1:N);
     [u,int_D,params]= ec_model2(N,Y(k,1:N),ecini,ts,sim_data,f,'plot_soln',0);
     %[x,y,dx,dy,u,int_D,params]=feval(scheme,N,ts,L,Y(k,1:N),Y_target,Omega,Mx,My,u_ext,f_src);
     for tc=1:length(ts)
      J(k,tc) =sum(sum(((u(:,:,tc)-u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2;
     end
     J_new(k)=dt*sum(J(k,1:length(ts)))+beta*int_D/2;
  end 
  RelError(k)=abs(J_new(k)-J_new(k-1))/abs(J_new(k));  
end
close all;
if ispc, ms = '\'; else ms = '/'; end
ad = cd;
cd(['..',ms,'Porous_Adjoint_Deterministic']);

for i=1:Mx
  for j=1:My
    D_target(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y_target(:));
    D(i,j)=feval('Dcoeff',sim_data.xcs(i),sim_data.ycs(j),N,Y(k,:));
 end
end 
  
figure;
plot(1:k,J_new(1:k), 'LineWidth',2);
title('Convergence of the cost functional J')
xlabel('Iterations');ylabel('Cost Functional');
print(gcf,[['Figures_Deterministic_N',num2str(N),num2str(log10(tol)),'tol'],ms,['fig1_',num2str(Mx)]],'-depsc')

figure;
plot(1:k,log10(J_new(1:k)), 'LineWidth',2);
title('Convergence of Log10(Cost Functional J)')
xlabel('Iterations');ylabel('Log10(Cost Functional)');
print(gcf,[['Figures_Deterministic_N',num2str(N),num2str(log10(tol)),'tol'],ms,['fig2_',num2str(Mx)]],'-depsc')

figure;
plot(1:k,Y(1:k,:));
title(['N= ',num2str(N),' trajectories of Y for a ',num2str(Mx),'X',num2str(My),' spatial mesh'])
xlabel('Iterations');ylabel('Y values');
print(gcf,[['Figures_Deterministic_N',num2str(N),num2str(log10(tol)),'tol'],ms,['fig3_',num2str(Mx)]],'-depsc')

figure;
plot(sim_data.xcs,u_target(:,My/2+1,end),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,u(:,My/2+1,end),'-b','LineWidth',2);hold off;
title('Crossections of target solution(red) versus estimated solution(blue) ')
print(gcf,[['Figures_Deterministic_N',num2str(N),num2str(log10(tol)),'tol'],ms,['fig4_',num2str(Mx)]],'-depsc')

figure;
plot(sim_data.xcs,D_target(:,My/2+1),'-r','LineWidth',2);hold on;
plot(sim_data.xcs,D(:,My/2+1),'-b','LineWidth',2);hold off;
title('Crossections of target diffusion(red) versus estimated diffusion(blue) ')
print(gcf,[['Figures_Deterministic_N',num2str(N),num2str(log10(tol)),'tol'],ms,['fig5_',num2str(Mx)]],'-depsc')

