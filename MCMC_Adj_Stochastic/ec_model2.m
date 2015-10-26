function [ecs,int_D,params] = ec_model2(N,param_vec,ecini,ets,g,f,varargin)

%  param_vec-vector of parameters, currently param_vec = [D,kp], just two
%    though it can be adjusted as need be so long as you adjust the program
%    accordingly, in particular the subprogram advance_ecs would need to be
%    changed.
%  ecini-initial epithelial cell density in the grid
%  ets-experimental times at which we want to know the epithelial cell 
%    density in the grid
%  g-grid info including
%    g.nx-number of cells in x-direction
%    g.ny-number of cells in y-direction
%    g.dx-size of cells in x-direction (we assume uniform spacing)
%    g.dy-size of cells in y-direction (again uniform spacing assumed)
%  Note:  The info inside g is later augmented with the appropriate time
%    step to use during each step in the numerical process (changes with
%    time).

%%  Extract info from param_vec
%  D-diffusion coefficient, corresponds to diffusion rate of cells
%  kp-growth or proliferation rate of epithelial cells
%  hc-hill coefficient, power in the hill term in the diffusion
%  Note:  To make the hill coefficient a free parameter, you need to swap
%  the commented and uncommented lines below for hc.  You would also need
%  to make the right changes in the MCMC program that called this program.
% params.D = param_vec;
[x_Dx,y_Dx] = meshgrid(g.xes,g.ycs);
[x_Dy,y_Dy] = meshgrid(g.xcs,g.yes);
[x_D,y_D]=meshgrid(g.xcs,g.ycs);
params.Dx= ones(g.ny,g.nx+1)+x_Dx.^2+y_Dx.^2;
params.Dy= ones(g.ny+1,g.nx)+x_Dy.^2+y_Dy.^2;
params.D=ones(g.ny,g.nx)+x_D.^2+y_D.^2;
%params.Dx= zeros(g.ny,g.nx+1);
%params.Dy= zeros(g.ny+1,g.nx);
for i=1:N 
    params.Dx=params.Dx+cos(i*pi*x_Dx/0.1).*cos(i*pi*y_Dx/0.07).*param_vec(i); 
    params.Dy=params.Dy+cos(i*pi*x_Dy/0.1).*cos(i*pi*y_Dy/0.07).*param_vec(i);   
    params.D=params.D+cos(i*pi*x_D/0.1).*cos(i*pi*y_D/0.07).*param_vec(i); 
end
params.Dx=params.Dx*(1e-3);
params.Dy=params.Dy*(1e-3);
params.D=params.D*(1e-3);
int_D=sum(params.D(:).*params.D(:).*g.dx.*g.dy);
%params.Dx=params.Dx/N+ones(g.ny,g.nx+1)+x_Dx.^2+y_Dx.^2;
%params.Dy=params.Dy/N+ones(g.ny+1,g.nx)+x_Dy.^2+y_Dy.^2;
 % params.kp = param_vec(2);
  % params.hc = param_vec(3);
 % params.hc = 2;

  %%  Consts for good numerical solving and other options
  ecm.safety_net =.5;
  ecm.plot_soln = 0;
  
  for vac = 1:2:length(varargin)
    ecm.(varargin{vac}) = varargin{vac+1};
  end

  %%  Get maximum dt allowable for stability
  %  This includes a safety buffer fraction of 0.9.
  %temp_quant = g.dx^2*g.dy^2/(max(params.Dx(:))*(g.dx^2+g.dy^2));
  g.dx;
  g.dy;
  max(max(params.Dx(:)),max(params.Dy(:)));
  temp_quant = (g.dx^2+g.dy^2)/(8*max(max(params.Dx(:)),max(params.Dy(:))));
  maxdt = ecm.safety_net*temp_quant;
 
  
  ecs = zeros([size(ecini),length(ets)]);
  ecs(:,:,1) = ecini;
  %figure;title(num2str(param_vec(1)))
  my_plot(ecm.plot_soln,ecini);

  for etc = 1:length(ets)-1

    exp_dt = ets(etc+1)-ets(etc);
    nsteps = ceil(exp_dt/maxdt);
    comp_dt = exp_dt/nsteps;
    g.dt = comp_dt;

    ecs(:,:,etc+1) = ecs(:,:,etc);
    for ctc = 1:nsteps
      ecs(:,:,etc+1) = advance_ecs(ecs(:,:,etc+1),g,params,f);
    end
    
    my_plot(ecm.plot_soln,ecs(:,:,etc+1));

  end
end

function my_plot(plot_soln,ecs)
 if plot_soln
   surf(ecs);
   zlim([0,1]);
   hold on; 
   contour(ecs,[0.5 0.5]);
   hold off;
   pause;
  end
end

function ecs = advance_ecs(ecs,g,params,f)
  
  grad_x = [(ecs(1,:)-0)./(g.dx/2);...
    diff(ecs,1,1)./g.dx;...
    (0-ecs(end,:))./(g.dx/2)];
  
  grad_y = [(ecs(:,1)-0)./(g.dy/2),...
    diff(ecs,1,2)./g.dy,...
    (0-ecs(:,end))./(g.dy/2)];
 
   ecs = ecs+g.dt*((diff(params.Dy.*grad_x,1,1)./g.dx+...
    diff(params.Dx.*grad_y,1,2)./g.dy)+f);

end
