function NEC_MCMC1(varargin)

  %  Performs MCMC

  %%  Plotting information (for later)
  clear all
  if matlabpool('size') == 0, matlabpool; end
  color = {'c-' 'b-' 'g-' 'r-' 'm-'};

  set(0,'DefaultAxesFontSize',20);
  set(0,'defaultaxeslinewidth',1);
  set(0, 'defaultlinelinewidth',   0.8);
  set(0, 'defaultpatchlinewidth',  0.7);

  %%  Initial guess and related information
  %  Number of independent variables
  N =5;
 
  %  Number of chains to use in parallel tempering
  Nc = 3;
  
  %  Measurement/likelihood standard deviation.  Note:  It will generally
  %  be the trend that if you lower mstd, you should also lower the
  %  baseline_jump size below.
  %mstd =1e-4; %for 2 vars
  %mstd =3e-5; %for 3 vars
  mstd =0.9e-5; %for 5 vars
  %  Baseline jump size.  The chain with the "lowest temperature" i.e. with
  %  the smallest jumps uses this baseline jump size.
  %baseline_jump = 0.001;
  %baseline_jump = 0.0001;
  
  %  The fractional jump sizes (jump_spreads) and actual jump sizes for the
  %  other chains are as follows.  Note that jump_spreads starts at 100%
  %  of the baseline value for the first chain and then moves upwards
  %  (sqrt(3)*100% next I believe)
  %jump_spreads = 1./(sqrt((1/3).^(([1:Nc]-1))));
  %jump_sizes = baseline_jump.*jump_spreads;
  
  %  To allow chains with larger jump sizes to actually explore space more,
  %  we have to adjust their corresponding likelihood functions to have
  %  wider effective standard deviations.  The following vector takes care
  %  of that.  Notice that likelihood_adj(1) = 1 always so the first
  %  chain's likelihood is never adjusted.
  %likelihood_adj = (1./jump_spreads);
  %beta=[1,0.34,0.1]  %for 5 vars
  beta=[1,0.345,0.1]
  jump_sizes=0.00155./sqrt(beta) %for 5 vars
  %beta=[1,0.2,0.05] %for 2 vars
  %jump_sizes=0.016./sqrt(beta)  %for 2 vars
  %beta=[1,0.3,0.1] %for 3 vars
  %jump_sizes=0.009./sqrt(beta)  %for 3 vars
  %  How many runs we make before checking whether or not we should swap
  %  chains.  If sims_per_swap = 100, we swap chains once every 100
  %  simulations.
  sims_per_swap = 50; %for 2,3,5 vars
%ranges for prior distribution [ 1/bnd to bnd ]
  bnd = ones(N,1);
  %  Initial guess for parameter values/parameter values used when using
  %  simulated experimental measurements.
  for i=1:N
    X_ini(1,i) = 0.55;
  end

  %  Limits on how far we explore with our guesses.  min_frac = 0.1 and
  %  max_frac = 10 means we keep our guesses within 10% and 1000% of X_ini.
  %  For instance if X_ini = [0.5,1,0.5,2,0.1], guesses will vary (think of
  %  this component-wise) between x_min = [0.05,0.1, 0.05,0.2,0.01] and x_max
  %  = [5, 10, 5, 20, 1].
  min_frac = 0.1;
  max_frac = 10;

  %  Number of guesses to make.  For every guess we try, we make a
  %  corresponding simulation (when interpolation is off).
  nruns =100000;

  %  First accepted guesses, just use X_ini for the first accepted guess
  %  for all chains (can change this later if we want).  repmat takes X_ini
  %  and changes it into [X_ini;X_ini;X_ini] when Nc = 3...that is it
  %  reproduces X_ini 3 times and lays it into a matrix 3 times in the row
  %  direction and 1 time in the column direction.
  x = repmat(X_ini,Nc,1);

  %  Find the first acccepted guess's corresponding likelihood value.  At the
  %  same time this initial call also reads in the experimental data
  %  (simulated or real) and creates an interpolation structure for use
  %  later.  IMPORTANT, this call and the call within the loop below must be
  %  consistent with each other.  E.g. if you set 'sim_data' to 1 here, you
  %  should also set 'sim_data' to 1 below, similarly with 'interpolate'
 % for cc = 1:Nc
 %   [fx(cc), interp_structure, exper_data] = likelihood_ecs1(N,x(cc,:),...
 %     [],[],'interpolate',0,'sim_data',1,'exp_plot',0,'plot_solns',0,...
 %     'meas_stddev',mstd,'params_baseline',X_ini');
 % end
 for cc = 1:Nc
    [fx(cc), interp_structure, exper_data] = likelihood_ecs1(N,x(cc,:),...
      [],[],'interpolate',0,'sim_data',1,'exp_plot',0,'plot_solns',0,...
      'meas_stddev',mstd);
 end
  accept = 0;
  accept_count = 0;
  xs = cell(1,Nc);
  fs = cell(1,Nc);
  nruns = sims_per_swap*ceil(nruns/sims_per_swap);
  for cc = 1:Nc
    xs{cc} = zeros(nruns,N);
    fs{cc} = zeros(nruns,1);
  end
  fprintf('Note:  Using %f simulations.\n',nruns);
  
  %  Initialize unaccepted guess variable:
  y = x;
  fy = fx;
  accept = zeros(1,Nc);
  swap_accepted = zeros(1,Nc-1);
  swap_attempted = zeros(1,Nc-1) ;
  accept_count = zeros(1,Nc);

  %  Outer for loop...we loop through sims_per_swap simulations
  for swapc=1:ceil(nruns/sims_per_swap)
    
    parfor cc = 1:Nc
     
      for simc = 1:sims_per_swap
        
        myc = (swapc-1)*sims_per_swap+simc;
        %  Use transition probability to propose a new guess as to suitable
        %  parameters for the model.  Note the new way I put in is faster.
       y(cc,:) = x(cc,:).*exp(normrnd(0,jump_sizes(cc)*ones(1,N)));
       %y(cc,:) = x(cc,:)+normrnd(0,jump_sizes(cc)*ones(1,N));
        %y(cc,:) = x(cc,:).*exp(randn(1,N)*jump_sizes(cc));
%         for j=1:N
%           y(j)=x(j)*exp(normrnd(0,.001)); %"x + jump"
%         end
       y_row= y(cc,:);     
            for i=1:N
                if y_row(i) > bnd(i)
                    y_row(i) = y_row(i)/bnd(i)^2;
                    %y_row(i) = y_row(i)-bnd(i);
                elseif y_row(i) < 1/bnd(i)
                   y_row(i) = y_row(i)*bnd(i)^2;
                end;
            end;
        y(cc,:)=y_row;    
        %  Assess how likely it is that the new guess for the parameters
        %  actually produced the (simulated or real) experimental measurements.
        %  Again note consistency of this call with what is above for
        %  'sim_data' and 'interpolate'
        [fy(cc),~,~] = likelihood_ecs1(N,y(cc,:),interp_structure,exper_data,...
          'interpolate',0,'sim_data',1,'exp_plot',0,'plot_solns',0,...
          'meas_stddev',mstd);
        
        %  Relative likelihood that the new guess actually produced the
        %  data compared to the previous accepted guess.  Note, the larger
        %  likelihood_adj is, the more likely it is that we'll accept the
        %  guess.
       %h = min(1, exp((fy(cc)-fx(cc)).*likelihood_adj(cc)));
       h = min(1, exp((fx(cc)-fy(cc)).*beta(cc)));
       % h = min(1, exp((fy(cc)-fx(cc))))
        %  Use U to {\it sometimes} accept new guesses that have lower
        %  likelihoods than the previously accepted guess
        U = rand;
        
        %  Either accept the new guess or keep the old guess by duplicating
        %  it/putting it into the markvo chain another time (by not doing
        %  anything to x by skipping the below if statement, we automatically
        %  copy the previous accepted guess into being the next accepted
        %  guess...see below)
        if U <= h
          x(cc,:) = y(cc,:);
          fx(cc) = fy(cc);
          accept(cc) = accept(cc) + 1;
        end
        accept_count(cc) = accept_count(cc) + 1;
        
        %  Put the most recent accepted guess into the markov chain as well as
        %  the corresponding likelihood
        xs{cc}(myc,:) = x(cc,:);  % Save chain into vector
        fs{cc}(myc) = fx(cc);
      end
    end
    
    for cc = Nc:-1:2
      if exp((beta(cc)-beta(cc-1))*(fx(cc)-fx(cc-1))) > rand(1)
        temp = x(cc,:);
        x(cc,:) = x(cc-1,:);
        x(cc-1,:) = temp;
        temp = fx(cc);
        fx(cc) = fx(cc-1);
        fx(cc-1) = temp;
        swap_accepted(cc-1)=swap_accepted(cc-1)+1;
      end
      swap_attempted(cc-1)=swap_attempted(cc-1)+1;
      %fprintf('%g/%g\n',swapc,ceil(nruns/sims_per_swap))
    end
      
  end
matlabpool close
  xs = xs{1}';
  fs = fs{1};

filename=strcat(num2str(N),'vars_',num2str(nruns/1000),'kruns_',num2str(sims_per_swap),'swaps.mat')
save(filename)
 
 
  %  Acceptance ratios and expected value.
 % acceptance = accept./accept_count
 % swap_rate=swap_accepted./swap_attempted
 % expected_value = mean(xs,2)
  
 % myinds = round(linspace(1,nruns,min(1000,nruns)));

 % figure;
 % plot(myinds, xs(:,myinds));

 % figure;
 % plot(myinds,fs{1}(myinds),'r',myinds,fs{2}(myinds),'b',myinds,fs{3}(myinds),'g');

 % figure;
 % plot(myinds,exp(fs{1}(myinds)),'r',myinds,exp(fs{2}(myinds)),'b',myinds,exp(fs{3}(myinds)),'g');
 
 % figure;
 % for xc = 1:N
 %   for xc2 = xc:N
 %     subplot(N,N,sub2ind([N,N],xc,xc2));
 %     plot(xs(xc,myinds),xs(xc2,myinds),'.');
 %     axis tight
 %   end
 % end

 % figure
 % for xc = 1:size(xs,1)
 %   subplot(1,size(xs,1),xc)
 %   hist(xs(xc,:))
 % end

end