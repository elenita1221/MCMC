  clear;
  vars={'xs','fs','accept','accept_count','swap_accepted','swap_attempted','nruns','N'};
  load('5vars_100kruns_50swaps.mat',vars{:});
  
  acceptance = accept./accept_count
  swap_rate=swap_accepted./(swap_attempted+swap_accepted)
  expected_value = mean(xs,2)'
  covariance_matrix=cov(xs')
  std_array=sqrt(diag(covariance_matrix))'
  
  myinds = round(linspace(1,nruns,min(1000,nruns)));
  close all;
  figure;
  %plot(myinds, xs(:,myinds));
  plot(myinds, xs(1,myinds),'r',myinds, xs(2,myinds),'g',myinds, xs(3,myinds),'b',...
      myinds, xs(4,myinds),'y',myinds, xs(5,myinds),'m');

 figure;
 plot(myinds,fs(myinds),'r');
 %plot(myinds,fs{1}(myinds),'r',myinds,fs{2}(myinds),'b',myinds,fs{3}(myinds),'g');

 figure;
 plot(myinds,exp(fs(myinds)),'r');
 %plot(myinds,exp(fs{1}(myinds)),'r',myinds,exp(fs{2}(myinds)),'b',myinds,exp(fs{3}(myinds)),'g');
 
  figure;
  for xc = 1:N
    for xc2 = xc:N
      subplot(N,N,sub2ind([N,N],xc,xc2));
      plot(xs(xc,myinds),xs(xc2,myinds),'.');
      set(gca,'FontSize',8);
      axis tight
    end
  end

  figure; 
  for xc = 1:size(xs,1)
    subplot(1,size(xs,1),xc)
    [n,X]=hist(xs(xc,:),20);
    bar(X,n/sum(n));
    set(gca,'XLim',[0.2 0.8],'FontSize',8);
  end  