function [ln_ml_val,interp_struc,exp_data] = ...
  likelihood_ecs1(N,param_vals,interp_struc,exp_data,varargin)

  %  param_vals-the parameter values at which we wish to find the 
  %  interp_struc-if you have no structure for interpolating, hand in []
  %    and the program will automatically find it for you
  %  exp_data-the corresponding experimental data, over all time
  %  largs-nice consts to have around
  if ispc, ms = '\'; else, ms = '/'; end
  largs.interpolate = 0;
  largs.min_frac = 0.1;
  largs.max_frac = 10;
  largs.sim_data = 1;
  largs.expdir = '';
  largs.nx = 10;
  largs.ny = 10;
  largs.ets = linspace(0,0.02,20);
 for i=1:N
  Y(i)=0.5;   
  largs.params_baseline(i)=Y(i);
 end 
  largs.meas_stddev = 1e-5;
  largs.exp_plot = 0;
  largs.plot_solns = 0;
  f=zeros(largs.nx,largs.ny);
  for vac = 1:2:length(varargin)
    largs.(varargin{vac}) = varargin{vac+1};
  end
 
  if isempty(exp_data)
    exp_data.g.nx = largs.nx;
    exp_data.g.ny = largs.ny;
    if largs.sim_data == 1
      exp_data.ts = largs.ets;
      exp_data.g.xl = -0.05;
      exp_data.g.xr = 0.05;
      exp_data.g.dx = (exp_data.g.xr-exp_data.g.xl)./exp_data.g.nx;
      exp_data.g.yl = -0.035;
      exp_data.g.yr = 0.035; 
      exp_data.g.dy = (exp_data.g.yr-exp_data.g.yl)./exp_data.g.ny;
      
      exp_data.g.xes = ...
        linspace(exp_data.g.xl,exp_data.g.xr,exp_data.g.nx+1);
      exp_data.g.yes = ...
        linspace(exp_data.g.yl,exp_data.g.yr,exp_data.g.ny+1);
      exp_data.g.xcs = (exp_data.g.xes(1:end-1)+exp_data.g.xes(2:end))./2;
      exp_data.g.ycs = (exp_data.g.yes(1:end-1)+exp_data.g.yes(2:end))./2;
      
      [xM,yM] = meshgrid(exp_data.g.xcs,exp_data.g.ycs);
     
      ecini = double((xM.^2+yM.^2)<0.01^2);
    % ecini = double((xM.^2+yM.^2)>0.01^2);
      
      exp_data.ecs = ec_model1(N,largs.params_baseline,ecini,exp_data.ts,...
        exp_data.g, f,'plot_soln',0);
      %exp_data.noise=ones([size(ecini),length(exp_data.ts)])*1e-4;
      %exp_data.ecs=exp_data.ecs+exp_data.noise;
    else
%       keyboard
      temp_data = load('../NEC Data/KalmanDataNEC_nx10_ny10.mat');
      exp_data.ts = temp_data.TimeExperiments;
      
      exp_data.g.nx = temp_data.XCells;
      exp_data.g.ny = temp_data.YCells;
      exp_data.g.xl = -temp_data.exp_ext(1)/2;
      exp_data.g.xr = temp_data.exp_ext(1)/2;
      exp_data.g.dx = (exp_data.g.xr-exp_data.g.xl)./exp_data.g.nx;
      exp_data.g.yl = -temp_data.exp_ext(2)/2;
      exp_data.g.yr = temp_data.exp_ext(2)/2;
      exp_data.g.dy = (exp_data.g.yr-exp_data.g.yl)./exp_data.g.ny;
      
      exp_data.g.xcs = squeeze(temp_data.ExpDataX(1,1,:));
      exp_data.g.ycs = squeeze(temp_data.ExpDataY(1,:,1));
      
      ecini = squeeze(temp_data.ExpData(1,:,:));
      exp_data.ecs = zeros([size(ecini),length(exp_data.ts)]);
      for tc = 1:length(exp_data.ts)
        exp_data.ecs(:,:,tc) = squeeze(temp_data.ExpData(1+(tc-1)*1,:,:));
        if largs.exp_plot == 1
          surf(exp_data.ecs(:,:,tc));
          hold on;
          contour(exp_data.ecs(:,:,tc),[0.5,0.5]);
          zlim([0,1]);
          hold off;
          pause
        end
      end
      
%       error('Real data not yet implemented!\n');
    end
  end
  
  if largs.interpolate == 0
    my_soln = ec_model1(N,param_vals,exp_data.ecs(:,:,1),exp_data.ts,exp_data.g,f,...
      'plot_soln',largs.plot_solns);
  else
    addpath(['..',ms,'interpolation_folder']);
    nout = exp_data.g.nx*exp_data.g.ny*length(exp_data.ts);
    if isempty(interp_struc)
      interp_struc = get_interpolant('myf', @vec_to_vec_func, ...
        'range', largs.params_baseline*[largs.min_frac,largs.max_frac],...
        'nout', nout, 'plot_things', 0, ...
        'MinDepth', 0, 'MaxDepth', 3, ...
        'cell_of_addl_args', ...
        {exp_data.ecs(:,:,1),exp_data.ts,exp_data.g,nout,'plot_soln',largs.plot_solns});
    end
    my_soln = zeros([size(exp_data.ecs),size(param_vals,2)]);
    for dc = 1:nout
      [i1,i2,i3] = ind2sub(size(exp_data.ecs),dc);
      interp_struc.selectOutput = dc;
      my_soln(i1,i2,i3,:) = spinterp(interp_struc, param_vals')';
    end
  end
  
  ln_ml_val = sum((my_soln(:)-exp_data.ecs(:)).^2)./2./...
   largs.meas_stddev^2./exp_data.g.nx./exp_data.g.ny./...
    length(exp_data.ts);
%   return
%   
%   global z;
%   if construct_interp == 1
%     addpath(['..',ms,'interpolation_folder']);
%     param_vals_baseline = varargin{1}';
%     param_names = varargin{2};
%     [trash,comp_meas,exp_meas,dirname,exp_data,norm_consts] = max_likelihood_prog(...
%       {mlpwi.mlp_arg_pairs{:},'sim_data_params',mlpwi.sim_data_params,...
%       'sim_ts',mlpwi.sim_ts},param_vals_baseline,param_names);
%     addpath('..\..\interpolation_folder');
%     z = get_interpolant('myf', @vec_to_vec_func, ...
%       'range', param_vals_baseline*[mlpwi.min_frac,mlpwi.max_frac],...
%       'nout', numel(comp_meas), 'plot_things', 0, ...
%       'MinDepth',0, 'MaxDepth', 1,...
%       'cell_of_addl_args',...
%       {{mlpwi.mlp_arg_pairs{:},'exp_data',exp_data},param_names,numel(comp_meas)});
% %     load z_struc;
%     z.exp_meas = exp_meas;
% %     temp = load([dirname{1},'Lung_data',ms,'Lung_run']);
%     z.nc_w_val = norm_consts{1};
%     fn = fieldnames(z.nc_w_val);
%     for fnc = 1:length(fn)
%       z.nc.(fn{fnc}) = z.nc_w_val.(fn{fnc}).val;
%     end
%     z.ml_mult = (exp_data.comp_domain.xe(2)-exp_data.comp_domain.xe(1)).*...
%       (exp_data.comp_domain.ye(2)-exp_data.comp_domain.ye(1))./...
%       ((exp_data.comp_domain.xe(end)-exp_data.comp_domain.xe(1))*...
%       (exp_data.comp_domain.ye(end)-exp_data.comp_domain.ye(1)));
%     z.nx = length(exp_data.comp_domain.xe);
%     z.ny = length(exp_data.comp_domain.ye);
% %     z.nz = temp.nz;
%     save('z_struc','z');
%   end
%   interp_soln = zeros([size(z.exp_meas),size(param_vals,2)]);
%   for dc = 1:z.d
%     [i1,i2,i3] = ind2sub(size(z.exp_meas),dc);
%     z.selectOutput = dc;
%     interp_soln(i1,i2,i3,:) = spinterp(z, param_vals')';
%   end
% %   parfor dc = 1:z.d
% %     ztemp = z;
% %     ztemp.selectOutput = dc;
% %     temp(dc,:) = spinterp(ztemp, param_vals')';
% %   end
% %   for dc = 1:z.d
% %     [i1,i2,i3] = ind2sub(size(z.exp_meas),dc);
% %     interp_soln(i1,i2,i3,:) = temp(dc);
% %   end
%   
%   for cc = 1:size(param_vals,2)
%     for ncc = 1:length(param_names)
%       z.nc.(param_names{ncc}) = param_vals(ncc,cc);
%       z.nc_w_val.(param_names{ncc}).val = param_vals(ncc,cc);
%     end
%     
%     delta = 0.2;
%     ml_constraint_sca_kca = 10000*(z.nc.sca/z.nc.kca-0.125).^2;
%     ml_constraint_cabar = 10000*...
%       (z.nc.cabar-sqrt((delta*((0.5/0.125)*z.nc.sca/z.nc.kca)^2-(z.nc.sca/z.nc.kca)^2)./(1-delta))).^2;
%     Rcah = 1/(1+((z.nc.sca/z.nc.kca)/z.nc.cabar)^2);
%     ml_constraint_healthy_steady_state = 100000*...
%       ((z.nc.km/z.nc.kmcp)<Rcah*z.nc.kcpml*z.nc.mrest/z.nc.kcp);
%     ml_constraint_n_kill_more = 100000*(z.nc.kn/z.nc.kncp>z.nc.nb*delta*Rcah*z.nc.kcpml/z.nc.kcp);
%     ml_constraint_n_m_decay_equally = 10000*(z.nc.km/z.nc.kmcp-z.nc.kn/z.nc.kncp).^2;
%     ml_constraint_cytokines_decay_equally = 10000*(z.nc.kcp/z.nc.kcpml-z.nc.kca/z.nc.kcaml)^2;
%     z.nc_w_val = get_kmb(z.nc_w_val,1,'Delta_cp',(0.5/0.125-1)*z.nc.sca/z.nc.kca,'cah',z.nc.sca/z.nc.kca);
%     ml_constraint_kmb = 100000*(z.nc.kmb-z.nc_w_val.kmb.val).^2;
%     
%     penalty = 100000;
%     
%     ml_val(cc) = z.ml_mult*sum(sum(sum((interp_soln(:,:,:,cc)-z.exp_meas).^2)));
%     temp2 = param_names;
%     for fc = 1:length(param_vals)
%       temp2{fc} = num2str(param_vals(fc,cc));
%     end
%     ml_val(cc) = ml_val(cc)+penalty*(max(param_vals(:,cc)<0))+...
%       ml_constraint_sca_kca+...
%       ml_constraint_cabar+...
%       ml_constraint_healthy_steady_state+...
%       ml_constraint_n_kill_more+...
%       ml_constraint_n_m_decay_equally+...
%       ml_constraint_cytokines_decay_equally+...
%       ml_constraint_kmb;
%     temp2 = {param_names{:};temp2{:}};
%     %   disp(temp2);
%   end
end

function varargout = vec_to_vec_func(varargin)
  [x_info{1:length(varargin)-1}] = varargin{1:end-1};
  xs = cell2mat(x_info);
  fprintf('# of func evals in this sparse grid generation: %g\n',...
    size(xs,1));
  cell_of_addl_args = varargin{end};
  ecini = cell_of_addl_args{1};
  ets = cell_of_addl_args{2};
  g = cell_of_addl_args{3};
  nout = cell_of_addl_args{4};
  varargin_for_ec_model = {cell_of_addl_args{5:end}};
  [out_cell{1:nout}] = deal(ones(size(xs,1),1));
  if matlabpool('size') == 0
    matlabpool 2;
  end
  comp_meas_mat = zeros(size(xs,1),nout);
  parfor xc = 1:size(xs,1)
    ecs = ec_model1(N,xs(xc,:),ecini,ets,g,f,...
      varargin_for_ec_model{:});
    fprintf('%g\n',xc);
    comp_meas_mat(xc,:) = ecs(:);
  end
  for nc = 1:nout
    out_cell{nc} = comp_meas_mat(:,nc);
  end
  varargout = out_cell;
end