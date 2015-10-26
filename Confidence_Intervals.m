clear all;
close all;
figure;
%line1=linspace(exp_val_MCMC(1)-std_dev_MCMC(1),exp_val_MCMC(1)+std_dev_MCMC(1),10000);
%plot(1000,line1,'r',1000,exp_val_MCMC(1),'xr')
%set(gca, 'XTick', [1000,100000])
%set(gca,'XLim',[0 100000]);
xlabel('Iterations');ylabel('Confidence Intervals');
%hold on

exp_val_SCKF=[0.5023, 0.4975, 0.5019, 0.4986, 0.5005];
std_dev_SCKF=[0.0945, 0.0966, 0.0842, 0.0790, 0.0754];

exp_val_EnKF=[0.5028, 0.5022, 0.4997, 0.4985, 0.4997];
std_dev_EnKF=[0.1038, 0.0982, 0.1008, 0.0885, 0.0801];

exp_val_KLSCKF=[0.5026, 0.4974, 0.5021, 0.4985, 0.5005];
std_dev_KLSCKF=[0.0956, 0.0970, 0.0843, 0.0792, 0.0755];

exp_val_KLEnKF=[0.5088, 0.4991, 0.5109, 0.4924, 0.5005];
std_dev_KLEnKF=[0.0500, 0.0500, 0.0500, 0.0500, 0.0500];

exp_val_MCMC=[0.4972,0.5057,0.4976,0.4984,0.5018];
std_dev_MCMC=[0.0400,0.0954,0.0916,0.0505,0.0161];

exp_val_Y1=[exp_val_SCKF(1),exp_val_EnKF(1),exp_val_KLSCKF(1),exp_val_KLEnKF(1),exp_val_MCMC(1)]
std_dev_Y1=[std_dev_SCKF(1),std_dev_EnKF(1),std_dev_KLSCKF(1),std_dev_KLEnKF(1),std_dev_MCMC(1)]

exp_val_Y2=[exp_val_SCKF(2),exp_val_EnKF(2),exp_val_KLSCKF(2),exp_val_KLEnKF(2),exp_val_MCMC(2)]
std_dev_Y2=[std_dev_SCKF(2),std_dev_EnKF(2),std_dev_KLSCKF(2),std_dev_KLEnKF(2),std_dev_MCMC(2)]

exp_val_Y3=[exp_val_SCKF(3),exp_val_EnKF(3),exp_val_KLSCKF(3),exp_val_KLEnKF(3),exp_val_MCMC(3)]
std_dev_Y3=[std_dev_SCKF(3),std_dev_EnKF(3),std_dev_KLSCKF(3),std_dev_KLEnKF(3),std_dev_MCMC(3)]

exp_val_Y4=[exp_val_SCKF(4),exp_val_EnKF(4),exp_val_KLSCKF(4),exp_val_KLEnKF(4),exp_val_MCMC(4)]
std_dev_Y4=[std_dev_SCKF(4),std_dev_EnKF(4),std_dev_KLSCKF(4),std_dev_KLEnKF(4),std_dev_MCMC(4)]

exp_val_Y5=[exp_val_SCKF(5),exp_val_EnKF(5),exp_val_KLSCKF(5),exp_val_KLEnKF(5),exp_val_MCMC(5)]
std_dev_Y5=[std_dev_SCKF(5),std_dev_EnKF(5),std_dev_KLSCKF(5),std_dev_KLEnKF(5),std_dev_MCMC(5)]

x=[1,2,3,4,5];
errorbar(x,exp_val_Y1,std_dev_Y1,'xr')
set(gca, 'XTick', [1,2,3,4,5]); set(gca,'XLim',[0 6]); set(gca,'YLim',[0.39 0.61]);
text(0.8, 0.45, 'SCKF', 'Color', 'k');
text(1.8, 0.45, 'EnKF', 'Color', 'k');
text(2.7, 0.45, 'KLSCKF', 'Color', 'k');
text(3.7, 0.45, 'KLEnKF', 'Color', 'k');
text(4.75, 0.45, 'MCMC', 'Color', 'k');
%title('Confidence Intervals for Y_1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
errorbar(x,exp_val_Y2,std_dev_Y2,'xb')
set(gca, 'XTick', [1,2,3,4,5]); set(gca,'XLim',[0 6]); set(gca,'YLim',[0.39 0.61]);
text(0.8, 0.45, 'SCKF', 'Color', 'k');
text(1.8, 0.45, 'EnKF', 'Color', 'k');
text(2.7, 0.45, 'KLSCKF', 'Color', 'k');
text(3.7, 0.45, 'KLEnKF', 'Color', 'k');
text(4.75, 0.45, 'MCMC', 'Color', 'k');
%title('Confidence Intervals for Y_2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
errorbar(x,exp_val_Y3,std_dev_Y3,'xg')
set(gca, 'XTick', [1,2,3,4,5]); set(gca,'XLim',[0 6]); set(gca,'YLim',[0.39 0.61]);
text(0.8, 0.45, 'SCKF', 'Color', 'k');
text(1.8, 0.45, 'EnKF', 'Color', 'k');
text(2.7, 0.45, 'KLSCKF', 'Color', 'k');
text(3.7, 0.45, 'KLEnKF', 'Color', 'k');
text(4.75, 0.45, 'MCMC', 'Color', 'k');
%title('Confidence Intervals for Y_3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
errorbar(x,exp_val_Y4,std_dev_Y4,'xm')
set(gca, 'XTick', [1,2,3,4,5]); set(gca,'XLim',[0 6]); set(gca,'YLim',[0.39 0.61]);
text(0.8, 0.45, 'SCKF', 'Color', 'k');
text(1.8, 0.45, 'EnKF', 'Color', 'k');
text(2.7, 0.45, 'KLSCKF', 'Color', 'k');
text(3.7, 0.45, 'KLEnKF', 'Color', 'k');
text(4.75, 0.45, 'MCMC', 'Color', 'k');
%title('Confidence Intervals for Y_4');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
errorbar(x,exp_val_Y5,std_dev_Y5,'xr')
set(gca, 'XTick', [1,2,3,4,5]); set(gca,'XLim',[0 6]); set(gca,'YLim',[0.39 0.61]);
text(0.8, 0.45, 'SCKF', 'Color', 'k');
text(1.8, 0.45, 'EnKF', 'Color', 'k');
text(2.7, 0.45, 'KLSCKF', 'Color', 'k');
text(3.7, 0.45, 'KLEnKF', 'Color', 'k');
text(4.75, 0.45, 'MCMC', 'Color', 'k');
%title('Confidence Intervals for Y_5');

