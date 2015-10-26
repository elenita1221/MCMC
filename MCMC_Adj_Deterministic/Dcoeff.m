function [D]=Dcoeff(x,y,N,param_vec)
D=1+x^1+y^2;
for k=1:N 
   D=D+cos(k*pi*x/0.1).*cos(k*pi*y/0.07).*param_vec(k);
end
D=D.*1e-3;   
