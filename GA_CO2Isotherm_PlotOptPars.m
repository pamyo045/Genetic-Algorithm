function fun = GA_CO2Isotherm_PlotOptPars
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
    
    % After multiple executions of this code, the best pars for each
    % pressure scenerio that I could extract (lowest SSR) are the
    % following:
    pars1={21.27,14.39,-4.77,-5.57,0.70,1.39,3.23,-1.26e+02};
    pars2={38.91,7.383,-4.24,-44.43,5.997,5.77,6.8,-233.2};
    pars3={10.46,4.955,0.4485,-5.433,10.26,0.217,24.06,-90.88};
    pars4={49.42,3.604,-3.698,-67.32,-5.171,5.85,10.58,-184.8};
    pars5={13.68,3.404,-1.895,-19.71,3.927,3.42,4.84,-180.2};
    
    for i=1:length(yKp_exp) 
        SSR = residual(yKp_exp{i});

        yModel = linspace(0,1,100);
        
        if i==1, optpars=pars1;
            elseif i==2, optpars=pars2;
            elseif i==3, optpars=pars3;
            elseif i==4, optpars=pars4;
            elseif i==5, optpars=pars5;
        end
        
        subplot(3,2,i);
        plot(yKp_exp{i}(:,1), yKp_exp{i}(:,2),'or',...
            yModel,Kp_model(yModel,optpars{:})),'-b';
        ylim([0 1.10*yKp_exp{i}(1,2)]);
        xlim([0 1]);
        xl=xlim;
        yl=ylim;
        legend('Experimental','Fitted Model');
        title(sprintf('System Pressure: %d atm',i));
        text(0.5*xl(2),0.5*yl(2),sprintf('SSR=%.4f',SSR(optpars)));
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
    
    function fit = fitness(args)
        fit = sumsqr(Kp_exp-Kp_mod(args{:}));
    end
    
    SSR=@fitness;
end

