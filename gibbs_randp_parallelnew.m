% Main file: where parameter estimation is implemented: Metropolis monte
% carlo with parallel tempering
% calls: fun_gibbs_new
%CHILD (reference ?)

function gibbs_randp_parallelnew(runnumber,foldername)
%e.g.gibbs_randp_parallelnew(1,'Outputtest')

global a %tdata y0 yind ydata stddata

%load datafile4.mat %only need this if you don't run using batch_gibbs.m

%number of chains
Nc = 3;
%number of parameters
Np = 4;
%number of param sets to save in one run
N =100000;

%distribution of Metropolis parameters
n = [1:Nc];
beta = (1/3).^(n-1);
epsilon = 0.3./sqrt(beta);
% To increase acceptance ratio: Decrease beta or epsilon, or increase stddata.


%baseline parameter values = A matrix values
pnorm = a;

%ranges for prior distribution [ 1/bnd to bnd ]
bnd = 10*ones(Np,1);


%##########################################################################
t = clock;
Njobs = 1;

%randomizing the generator
% rand('state',sum(100*clock));
% randn('state',sum(100*clock));

accepted = ones(1,Nc);
rejected = ones(1,Nc);
swap = zeros(1,Nc-1);
swattempt = 0;

warning off

pvec = single(zeros(Nc,Np,N));
Evec = single(zeros(Nc,N));

%infovec = single(zeros(14,N));
%gibbsdir = 'gibbs_output_widenew1/';


runnuminp = max(0,runnumber-Njobs);

%when runnumber==1
if runnuminp == 0
    p = ones(Nc,Np);
    
    %when runnumber > 1
else
    %load starting parameter values from the previous run
    gibbs_last_filename = [foldername '/gibbs_output_',num2str(runnuminp),'.mat'];
    
    load(gibbs_last_filename);
    p = pvec(:,:,end);
    p = double(p);    %%%%%%%%%%Not sure what this does
    
end

%info = zeros(1,14);

%randomizing initial condition    %%%%%Not sure what 'randomizing' means here
for kk = 1:Nc
    %slow version: use for defective case
    %  Ea(kk,1) = fun_gibbs_new(p(kk,:).*pnorm);  %x is input here
    
    %exact solution version (good for all cases except defective)
    Ea(kk,1) = exact_fun_gibbs_new(p(kk,:).*pnorm);  %x is input here (for exact_fun_gibbs_new)
    
end
j = 0;

%keep track of indeces where swaps occur
swapvector=zeros(2,N);

while j < N
    
    %run chain
    for ll = 1:Np  %run each chain Np times before swapping
        for kk = 1:Nc
            
            %%%%%This always outputs a positive number since p(kk) is initially 1
            pn = p(kk,:).*exp(randn(1,Np)*epsilon(kk));    %generate new set of parameters
            %perturbation is gaussian in the log space
            
            %randn(1,Np) 1xNp matrix of pseudo random values from standard normal
            %distribution (mu=0, sigma=1) (e.g. -3<entries<3)
            
            for i=1:Np
                if pn(i) > bnd(i)
                    pn(i) = pn(i)/bnd(i)^2;
                elseif pn(i) < 1/bnd(i)
                    pn(i) = pn(i)*bnd(i)^2;
                end;
            end;
            
            %slow version: use for defective case
            %Eb  = fun_gibbs_new(pn.*pnorm);        %calculate E(pn) %%%Eb=logPsum
            
            %exact solution approach (good for all cases except defective)
            Eb  = exact_fun_gibbs_new(pn.*pnorm);
            
            
            if exp(beta(kk)*(Ea(kk)-Eb))>rand(1)
                p(kk,:) = pn;
                Ea(kk) = Eb;
                accepted(kk) = accepted(kk) + 1;
            else
                rejected(kk) = rejected(kk) + 1;
            end;
        end
        j = j+1;
        pvec(:,:,j) = single(p);
        Evec(:,j) = single(Ea);
    end
    
    %try to swap chains
    swattempt = swattempt + 1;
    for kk = Nc:-1:2
        if exp((beta(kk)-beta(kk-1))*(Ea(kk)-Ea(kk-1))) > rand(1)
            Eb = Ea(kk);
            Ea(kk) = Ea(kk-1);
            Ea(kk-1) = Eb;
            pn = p(kk,:);
            p(kk,:) = p(kk-1,:);
            p(kk-1,:) = pn;
            swap(kk-1) = swap(kk-1)+1;
            swapvector(4-kk,j)=1;
        end
    end
    
    if mod(j,5000) == 0
        disp(sprintf('%d/%d',j,N));
        accratio = accepted./(accepted+rejected);
        swapratio = swap./(swap+swattempt);
    end
    
end

accratio = accepted./(accepted+rejected);
swapratio = swap./(swap+swattempt);

%save data in a file
gibbs_output_filename = [foldername '/gibbs_output_',num2str(runnumber),'.mat'];
save(gibbs_output_filename,'pvec','Evec','accratio','swapratio', 'N', 'swapvector');

disp(sprintf('Elapsed time:%f',etime(clock,t)))