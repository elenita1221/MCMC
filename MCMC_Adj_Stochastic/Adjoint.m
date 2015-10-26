clear all %Adjoint Parabolic Code for the Stochastic Case
N=5; nruns=50; eps=50000; beta=1e-6; tol=1e-7; eps=2*eps/3; % use 10 & 50 runs
Mx=10; My=10; %%Omega=[-0.05, 0.05, -0.035, 0.035];grid size
mySeed = 10; rng(mySeed); % only 20x20 grid
%k=1;Y(k,1:N) =0.6+randn(1,N)/100;  %initial guess for 1st iteration
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
E_u_target=zeros(Mx,My,length(ts));E_k_square=zeros(Mx,My);
E_u_target_square=zeros(Mx,My,length(ts));
noise=1e-3; %use noise 1e-1, 1e-3
for run_counter=1:nruns 
   Y_target(1:N,run_counter) = 0.5+noise.*randn(1,N);
   [u,~,~]= ec_model2(N,Y_target(:,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
   u_target(1:Mx,1:My,1:length(ts),run_counter)=u;
   E_u_target=E_u_target+u_target(:,:,:,run_counter);
   E_u_target_square=E_u_target_square+u_target(:,:,:,run_counter).*u_target(:,:,:,run_counter);
end
E_u_target=E_u_target/nruns; E_u_target_square=E_u_target_square/nruns;

for i=1:Mx
 for j=1:My
    cos_prod(i,j,1:N)=cos((1:N)*pi*sim_data.xcs(i)).*cos((1:N)*pi*sim_data.ycs(j))./(1e+3);
    cos_sum(i,j)=sum(cos_prod(i,j,:))./(1e+3);
 end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k3=1; sol3=cell(1,nruns); RelError3(1)=1e1; %J3 cost functional
E_u3=zeros(Mx,My,length(ts)); E_u3_square=zeros(Mx,My,length(ts));
for run_counter=1:nruns 
 Y_ini(1:N,run_counter) = 0.7+randn(1,N)*1e-2; 
 Y3(k3,1:N,run_counter)=Y_ini(1:N,run_counter);
 [u3,int_D,params]= ec_model2(N,Y3(k3,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
 sol3{run_counter}=u3;
for tc=1:length(ts)
  J(run_counter,k3,tc) =sum(sum(((u3(:,:,tc)-u_target(:,:,tc,run_counter)).^2).*sim_data.dx.*sim_data.dy))/2; %cost functional
end
J_new(run_counter,k3)=dt*sum(J(run_counter,k3,1:length(ts)))+beta*int_D/2;
E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3;
end
J3(k3)=mean(J_new(1:nruns,k3));
E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
quant3(k3)=dt*sum(sum(sum(((E_u3(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));

while RelError3(k3)>tol
  eps=3*eps/2; 
  k3=k3+1;
  E_u3=zeros(Mx,My,length(ts)); E_u3_square=zeros(Mx,My,length(ts));
 for run_counter=1:nruns 
  u3=sol3{run_counter};   
  [eta3,~,~]= ec_model2_adjoint(N,Y3(k3-1,1:N,run_counter),ecini,ts,sim_data,u3-u_target(:,:,:,run_counter),'plot_soln',0);
  for tc=1:length(ts)
     gradx_u = [diff(u3(:,:,tc),1,1)./sim_data.dx];
     grady_u = [diff(u3(:,:,tc),1,2)./sim_data.dy];
     gradx_eta=[diff(eta3(:,:,tc),1,1)./sim_data.dx];
     grady_eta=[diff(eta3(:,:,tc),1,2)./sim_data.dy];
     matrix1=gradx_u(:,:).*gradx_eta(:,:);
     matrix2=grady_u(:,:).*grady_eta(:,:);
     temp1=-sum(sum(cos_sum(1:Mx-1,:).*matrix1(:,:)))*sim_data.dx*sim_data.dy;
     temp2=-sum(sum(cos_sum(:,1:My-1).*matrix2(:,:)))*sim_data.dx*sim_data.dy;
     dJdY(k3-1,tc,run_counter)=dt*(temp1+temp2);
  end 
  for i=1:N
   dJdY_new(k3-1,i,run_counter)=sum(dJdY(k3-1,1:length(ts),run_counter))+beta*sum(sum(params.D(:,:).*cos_prod(:,:,i))).*sim_data.dx.*sim_data.dy;
  end
  Y3(k3,1:N,run_counter)=Y3(k3-1,1:N,run_counter)-eps*dJdY_new(k3-1,1:N,run_counter);
  [u3,int_D,params]= ec_model2(N,Y3(k3,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
  sol3{run_counter}=u3;
  for tc=1:length(ts)
    J(run_counter,k3,tc) =sum(sum(((u3(:,:,tc)-u_target(:,:,tc,run_counter)).^2).*sim_data.dx.*sim_data.dy))/2;
  end
  J_new(run_counter,k3)=dt*sum(J(run_counter,k3,1:length(ts)))+beta*int_D/2;
  E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3;
 end
 J3(k3)=mean(J_new(1:nruns,k3));
 E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
 quant3(k3)=dt*sum(sum(sum(((E_u3(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));
 
  while J3(k3)>=J3(k3-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u3=zeros(Mx,My,length(ts)); E_u3_square=zeros(Mx,My,length(ts));
    for run_counter=1:nruns
     Y3(k3,1:N,run_counter)=Y3(k3-1,1:N,run_counter)-eps*dJdY_new(k3-1,1:N,run_counter);
     [u3,int_D,params]= ec_model2(N,Y3(k3,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
     sol3{run_counter}=u3;
     for tc=1:length(ts)
      J(run_counter,k3,tc) =sum(sum(((u3(:,:,tc)-u_target(:,:,tc,run_counter)).^2).*sim_data.dx.*sim_data.dy))/2;
     end
     J_new(run_counter,k3)=dt*sum(J(run_counter,k3,1:length(ts)))+beta*int_D/2;
     E_u3=E_u3+u3; E_u3_square=E_u3_square+u3.*u3; 
    end 
    J3(k3)=mean(J_new(1:nruns,k3));
    E_u3=E_u3/nruns; E_u3_square=E_u3_square/nruns;
    quant3(k3)=dt*sum(sum(sum(((E_u3(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));
  end 
  
  RelError3(k3)=abs(J3(k3)-J3(k3-1))/abs(J3(k3));  
end
disp('J3 finished!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dJdY dJdY_new 
eps=50000; eps=2*eps/3; RelError4(1)=1e1;
k4=1; sol4=cell(1,nruns); c=0;  %J4 cost functional
E_u4=zeros(Mx,My,length(ts)); E_u4_square=zeros(Mx,My,length(ts));
for run_counter=1:nruns 
 Y4(k4,1:N,run_counter)=Y_ini(1:N,run_counter);
 [u4,int_D,params]= ec_model2(N,Y4(k4,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
 sol4{run_counter}=u4;
 E_u4=E_u4+u4; E_u4_square=E_u4_square+u4.*u4; E_k_square=E_k_square+params.D.*params.D;
end
E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
E_u4_square=E_u4_square/nruns;
for tc=1:length(ts)
  J4(k4,tc)=sum(sum(((E_u4(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
    +c*sum(sum(((E_u4_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
end  
J4_new(k4)=dt*sum(J4(k4,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
quant4(k4)=dt*sum(sum(sum(((E_u4(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));

old_E_u4=E_u4; old_E_u4_square=E_u4_square;
while RelError4(k4)>tol
     eps=3*eps/2;
     k4=k4+1;
     E_u4=zeros(Mx,My,length(ts)); E_k_square=zeros(Mx,My);
     E_u4_square=zeros(Mx,My,length(ts));
 for run_counter=1:nruns 
   u4=sol4{run_counter};
   [eta4,~,~]= ec_model2_adjoint(N,Y4(k4-1,1:N,run_counter),ecini,ts,sim_data,old_E_u4-E_u_target+c*u4.*(old_E_u4_square-E_u_target_square),'plot_soln',0);
   for tc=1:length(ts)
     gradx_u = [diff(u4(:,:,tc),1,1)./sim_data.dx];
     grady_u = [diff(u4(:,:,tc),1,2)./sim_data.dy];
     gradx_eta=[diff(eta4(:,:,tc),1,1)./sim_data.dx];
     grady_eta=[diff(eta4(:,:,tc),1,2)./sim_data.dy];
     matrix1=gradx_u(:,:).*gradx_eta(:,:);
     matrix2=grady_u(:,:).*grady_eta(:,:);
     temp1=-sum(sum(cos_sum(1:Mx-1,:).*matrix1(:,:)))*sim_data.dx*sim_data.dy;
     temp2=-sum(sum(cos_sum(:,1:My-1).*matrix2(:,:)))*sim_data.dx*sim_data.dy;
     dJdY(k4-1,tc,run_counter)=dt*(temp1+temp2);
   end 
   for i=1:N
    dJdY_new(k4-1,i,run_counter)=sum(dJdY(k4-1,1:length(ts),run_counter))+beta*sum(sum(params.D(:,:).*cos_prod(:,:,i))).*sim_data.dx.*sim_data.dy;
   end
   Y4(k4,1:N,run_counter)=Y4(k4-1,1:N,run_counter)-eps*dJdY_new(k4-1,1:N,run_counter);
   [u4,int_D,params]= ec_model2(N,Y4(k4,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
   sol4{run_counter}=u4;
   E_u4=E_u4+u4; E_k_square=E_k_square+params.D.*params.D;
   E_u4_square=E_u4_square+u4.*u4; 
 end
 E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
 E_u4_square=E_u4_square/nruns;
 for tc=1:length(ts)
  J4(k4,tc)=sum(sum(((E_u4(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
    +c*sum(sum(((E_u4_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
end  
J4_new(k4)=dt*sum(J4(k4,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
quant4(k4)=dt*sum(sum(sum(((E_u4(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));

  while J4_new(k4)>=J4_new(k4-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u4=zeros(Mx,My,length(ts)); E_k_square=zeros(Mx,My);
    E_u4_square=zeros(Mx,My,length(ts));
    for run_counter=1:nruns
     Y4(k4,1:N,run_counter)=Y4(k4-1,1:N,run_counter)-eps*dJdY_new(k4-1,1:N,run_counter);
     [u4,int_D,params]= ec_model2(N,Y4(k4,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
     %[x,y,dx,dy,u4,int_k,params]=feval(scheme,N,ts,L,Y4(k4,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol4{run_counter}=u4;
     E_u4=E_u4+u4; E_k_square=E_k_square+params.D.*params.D;
     E_u4_square=E_u4_square+u4.*u4; 
    end;
    E_u4=E_u4/nruns; E_k_square=E_k_square/nruns;
    E_u4_square=E_u4_square/nruns;
    for tc=1:length(ts)
     J4(k4,tc)=sum(sum(((E_u4(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
    +c*sum(sum(((E_u4_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
    end  
    J4_new(k4)=dt*sum(J4(k4,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
    quant4(k4)=dt*sum(sum(sum(((E_u4(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));  
  end
old_E_u4=E_u4; old_E_u4_square=E_u4_square;
RelError4(k4)=abs(J4_new(k4)-J4_new(k4-1))/abs(J4_new(k4));
end
disp('J4 finished!')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear dJdY dJdY_new E_k_square
eps=50000; eps=2*eps/3; RelError5(1)=1e1;
k5=1; sol5=cell(1,nruns); c=1;  %J5 cost functional
E_u5=zeros(Mx,My,length(ts)); E_u5_square=zeros(Mx,My,length(ts)); E_k_square=zeros(Mx,My);
for run_counter=1:nruns 
Y5(k5,1:N,run_counter)=Y_ini(1:N,run_counter);
[u5,int_D,params]= ec_model2(N,Y5(k5,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
%[x,y,dx,dy,u5,int_k,params]=feval(scheme,N,ts,L,Y5(k5,1:N,run_counter),Y_target(1:N,run_counter),Omega,Mx,My,u_ext,f_src);%SE parabolic solver
sol5{run_counter}=u5;
E_u5=E_u5+u5; E_u5_square=E_u5_square+u5.*u5; E_k_square=E_k_square+params.D.*params.D;
end;
E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
E_u5_square=E_u5_square/nruns;
for tc=1:length(ts)
  J5(k5,tc)=sum(sum(((E_u5(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
    +c*sum(sum(((E_u5_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
end  
J5_new(k5)=dt*sum(J5(k5,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
quant5(k5)=dt*sum(sum(sum(((E_u5(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));

old_E_u5=E_u5; old_E_u5_square=E_u5_square;
while RelError5(k5)>tol
     eps=3*eps/2;
     k5=k5+1;
     E_u5=zeros(Mx,My,length(ts)); E_k_square=zeros(Mx,My);
     E_u5_square=zeros(Mx,My,length(ts));
 for run_counter=1:nruns 
   u5=sol5{run_counter};
   [eta5,~,~]= ec_model2_adjoint(N,Y5(k5-1,1:N,run_counter),ecini,ts,sim_data,old_E_u5-E_u_target+c*u5.*(old_E_u5_square-E_u_target_square),'plot_soln',0);
   %[~,~,~,~,eta4,~,~]=feval('parabolic5pt_adjoint',N,ts,L,Y4(k4-1,1:N,run_counter),u4,u_target,Omega,Mx,My,old_E_u4-E_u_target+c*u4.*(old_E_u4_square-E_u_target_square));%AE parabolic solver
   for tc=1:length(ts)
     gradx_u = [diff(u5(:,:,tc),1,1)./sim_data.dx];
     grady_u = [diff(u5(:,:,tc),1,2)./sim_data.dy];
     gradx_eta=[diff(eta5(:,:,tc),1,1)./sim_data.dx];
     grady_eta=[diff(eta5(:,:,tc),1,2)./sim_data.dy];
     matrix1=gradx_u(:,:).*gradx_eta(:,:);
     matrix2=grady_u(:,:).*grady_eta(:,:);
     temp1=-sum(sum(cos_sum(1:Mx-1,:).*matrix1(:,:)))*sim_data.dx*sim_data.dy;
     temp2=-sum(sum(cos_sum(:,1:My-1).*matrix2(:,:)))*sim_data.dx*sim_data.dy;
     dJdY(k5-1,tc,run_counter)=dt*(temp1+temp2);
   end 
   for i=1:N
    dJdY_new(k5-1,i,run_counter)=sum(dJdY(k5-1,1:length(ts),run_counter))+beta*sum(sum(params.D(:,:).*cos_prod(:,:,i))).*sim_data.dx.*sim_data.dy;
   end
   Y5(k5,1:N,run_counter)=Y5(k5-1,1:N,run_counter)-eps*dJdY_new(k5-1,1:N,run_counter);
   [u5,int_D,params]= ec_model2(N,Y5(k5,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
   sol5{run_counter}=u5;
   E_u5=E_u5+u5; E_k_square=E_k_square+params.D.*params.D;
   E_u5_square=E_u5_square+u5.*u5; 
 end
 E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
 E_u5_square=E_u5_square/nruns;
 for tc=1:length(ts)
  J5(k5,tc)=sum(sum(((E_u5(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
    +c*sum(sum(((E_u5_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
 end  
 J5_new(k5)=dt*sum(J5(k5,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
 quant5(k5)=dt*sum(sum(sum(((E_u5(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy)));
 
 while J5_new(k5)>=J5_new(k5-1)
    eps=eps/10
    if eps<1e-15 disp('Algorithm stagnated');break; end;
    E_u5=zeros(Mx,My,length(ts)); E_k_square=zeros(Mx,My);
    E_u5_square=zeros(Mx,My,length(ts));
    for run_counter=1:nruns
     Y5(k5,1:N,run_counter)=Y5(k5-1,1:N,run_counter)-eps*dJdY_new(k5-1,1:N,run_counter);
     [u5,int_D,params]= ec_model2(N,Y5(k5,1:N,run_counter),ecini,ts,sim_data,f,'plot_soln',0);
     %[x,y,dx,dy,u5,int_k,params]=feval(scheme,N,ts,L,Y5(k5,1:N,run_counter),Y_target(:,run_counter),Omega,Mx,My,u_ext,f_src);
     sol5{run_counter}=u5;
     E_u5=E_u5+u5; E_k_square=E_k_square+params.D.*params.D;
     E_u5_square=E_u5_square+u5.*u5; 
    end;
    E_u5=E_u5/nruns; E_k_square=E_k_square/nruns;
    E_u5_square=E_u5_square/nruns;
    for tc=1:length(ts)
      J5(k5,tc)=sum(sum(((E_u5(:,:,tc)-E_u_target(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/2+...
     +c*sum(sum(((E_u5_square(:,:,tc)-E_u_target_square(:,:,tc)).^2).*sim_data.dx.*sim_data.dy))/4;
    end
    J5_new(k5)=dt*sum(J5(k5,1:length(ts)))+beta*sum(E_k_square(:).*sim_data.dx.*sim_data.dy)/2;%cost functional
    quant5(k5)=dt*sum(sum(sum(((E_u5(:,:,:)-E_u_target(:,:,:)).^2).*sim_data.dx.*sim_data.dy))); 
 end
old_E_u5=E_u5; old_E_u5_square=E_u5_square;
RelError5(k5)=abs(J5_new(k5)-J5_new(k5-1))/abs(J5_new(k5));
end
 
filename=strcat('All_J_N',num2str(N),'_',num2str(nruns),'runs',num2str(log10(tol)),'tol_',num2str(Mx),'x',num2str(My),'_noise',num2str(log10(noise)),'.mat')
save(filename)
