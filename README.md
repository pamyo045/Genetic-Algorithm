# Genetic-Algorithm
This code uses MATLAB's inherent function 'ga' to perform a genetic algorithm in order to find the minimum of a function - in this case the smallest sum of sqaured residual (SSR) between model parameters and experimental values within an excel file.

**The MATLAB code for GA_CO2Isotherm:**
```
function fun = GA_CO2Isotherm
    % XLSname: Binary - N2-CO2 with HISIV3000 silicalite combo graph.xlsx
    % XLSsheet: Kp exp - curve fits
    % XLSrange: A4:F22
    XLSname = 'Binary - N2-CO2 with HISIV3000 silicalite combo graph.xlsx';
    XLSsheet = 'Kp exp - curve fits';
    XLSrange = 'A4:F22';
    
    % extract experimental Kp values from Excel file and assign each col to
    % appropriate pressure trials
    XLSdata = readtable(XLSname,'sheet',XLSsheet,'range',XLSrange,'readvariablename',false);
    
    % organize data, variable index 1 of table = col 1 of spreadsheet
    yCO2_exp = XLSdata.(1); % store entire yCO2 col
    temp = XLSdata.(2); % assign temp array to col 2 data (i.e. 1atm data)
    yKp_exp_1atm = [yCO2_exp(~isnan(temp)), temp(~isnan(temp))]; % assign array without the zero Kp values along with corresponding yCO2 in a 2D array
    temp = XLSdata.(3);
    yKp_exp_2atm = [yCO2_exp(~isnan(temp)), temp(~isnan(temp))];
    temp = XLSdata.(4);
    yKp_exp_3atm = [yCO2_exp(~isnan(temp)), temp(~isnan(temp))];
    temp = XLSdata.(5);
    yKp_exp_4atm = [yCO2_exp(~isnan(temp)), temp(~isnan(temp))];
    temp = XLSdata.(6);
    yKp_exp_5atm = [yCO2_exp(~isnan(temp)), temp(~isnan(temp))];
    
    yKp_exp = {yKp_exp_1atm yKp_exp_2atm yKp_exp_3atm yKp_exp_4atm yKp_exp_5atm};
    figure;
    
    for i=1:length(yKp_exp) 
        SSR = residual(yKp_exp{i});

        options = gaoptimset('vectorized','off','TolFun',1e-9);
        [fun,fval] = ga(SSR,8,[],[],[],[],[],[],nlconsy1(yKp_exp{i}),options);
        fprintf('Trial: %d atm \nSSR=%.4g \n[B1,B2,B3,C1,C2,C3,gamma,lambda]=[%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g,%.4g]\n',...
            i,fval,fun(1),fun(2),fun(3),fun(4),fun(5),fun(6),fun(7),fun(8));

        for j=1:length(fun)
            optpars{j}=fun(j);
        end

        yModel = linspace(0,1,100);
        
        subplot(3,2,i);
        plot(yKp_exp{i}(:,1), yKp_exp{i}(:,2),'or',...
            yModel,Kp_model(yModel,optpars{:})),'-b';
        ylim([0 1.10*yKp_exp{i}(1,2)]);
        xlim([0 1]);
        xl=xlim;
        yl=ylim;
        legend('Experimental','Fitted Model');
        title(sprintf('System Pressure: %d atm',i));
        text(0.5*xl(2),0.5*yl(2),sprintf('SSR=%.4f',fval));
    end
end

function fun = Kp_model(y1,b1,b2,b3,c1,c2,c3,gamma,theta)
    fun=((1-y1).*(b1./abs(y1+gamma)+...
        b2.*exp(theta.*y1)+b3)+...
        y1.*(c1./abs(y1+gamma)+c2.*exp(theta.*y1)+c3));
end

function SSR = residual(yKp_exp)
    y1=yKp_exp(:,1);
    Kp_exp=yKp_exp(:,2);

    Kp_mod = @(b1,b2,b3,c1,c2,c3,gamma,theta)Kp_model(y1,b1,b2,b3,c1,c2,c3,gamma,theta);
    
    function fit = fitness(x)
        b1=x(:,1); b2=x(:,2); b3=x(:,3);
        c1=x(:,4); c2=x(:,5); c3=x(:,6);
        gamma=x(:,7); theta=x(:,8);
        args={b1,b2,b3,c1,c2,c3,gamma,theta};
        
        fit = sumsqr(Kp_exp-Kp_mod(args{:}));
    end
    
    SSR=@fitness;
end

function fun = nlconsy1(yKp_exp)
    y1=yKp_exp(:,1);
    Kp_exp=yKp_exp(:,2);
    
    Kp_mod = @(b1,b2,b3,c1,c2,c3,gamma,theta)Kp_model(y1,b1,b2,b3,c1,c2,c3,gamma,theta);
    
    function [c, ceq] = nlconsx(x)
        b1=x(:,1); b2=x(:,2); b3=x(:,3);
        c1=x(:,4); c2=x(:,5); c3=x(:,6);
        gamma=x(:,7); theta=x(:,8);
        args={b1,b2,b3,c1,c2,c3,gamma,theta};
        
        c = [sumsqr(Kp_exp-Kp_mod(args{:}))-0.2];
        ceq = [];
    end
    
    fun=@nlconsx;
    
end
```

**Example of a MATLAB Command Window Session:**
![](https://github.com/pamyo045/Genetic-Algorithm/blob/master/Image1.png)

**Example of the MATLAB code output:**
![](https://github.com/pamyo045/Genetic-Algorithm/blob/master/Image2.png)
