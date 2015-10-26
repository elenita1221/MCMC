%This file evaluates the objective function

function logPsum  = exact_fun_gibbs_new(x)

global tdata y0 ydata stddata yind

try
    %exact solution to the system with parameters x
    [y] = integrate_soln_mathematica(tdata,y0,x);
    
    ny = length(find(yind));
    nt = length(tdata); 
    
    logP = zeros(nt,ny);
    
    
    if yind(1)==1  %%this will be the case for yind=[1,1] or [1,0]
        
        for i = 1:nt % nt=# data points
            
            %computes difference between data and trajectory
            
            for j = find(yind)     %1:ny

         %%%%%% why isn't this  (2*stddata^2) ? (see below)
          logP(i,j) = (y(i,j)-ydata(i,j))^2/2/stddata(i,j)^2;
          
                %%%%%% why isn't this  (2*stddata^2) ?
            end
        end
        
    else %for the case yind=[0,1]
        
        for i = 1:nt % nt=# data points
            
            %want the second column of y and the 1st (only) column of ydata
            %and stdata
            logP(i,1) = (y(i,2)-ydata(i,1))^2/2/stddata(i,1)^2;
            
        end
    end
    
    %short version
% logP = (y(index,:)-ydata).^2./2./stddata.^2;
    
    logPsum = sum(sum(logP));
    
catch
    disp('Integration failed. Assigning infinite energy.');
    %keyboard
    logPsum = Inf;
end


%function that finds the exact solution (from mathematica) and evaluates it at times ts (tdata) 
    function ys = integrate_soln_mathematica(ts,y0,a)
        a11 = a(1); a12 = a(2); a21 = a(3); a22 = a(4);
        y10 = y0(1); y20 = y0(2);
        
        y1s = (-a11*exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + a22*exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)*exp(1/2*(a11 + a22 ...
            - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + a11*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            - a22*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)*exp(1/2*(a11 + a22+ sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            - 2*a12 *exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*...
            y20 + 2*a12*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y20)/(2*sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2));
        
        y2s = (-2*a21*exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + 2*a21*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)) .*ts)*y10...
            + a11*exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)*y20...
            - a22*exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)*y20...
            + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)* exp(1/2*(a11 + a22 - sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)*y20...
            - a11*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)*y20...
            + a22*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)*y20...
            + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)*exp(1/2*(a11 + a22 + sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2)).*ts)...
            *y20)/(2*sqrt(a11^2 + 4*a12*a21 - 2*a11*a22 + a22^2));
        ys = [y1s',y2s'];
    end
end
