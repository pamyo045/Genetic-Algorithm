function fits = ps_fitting(varargin)
    % xlsName: Binary - N2-CO2 with HISIV3000 silicalite combo graph.xlsx
    % xlsSheet: Kp exp - curve fits
    if ~isempty(varargin)
        initialPars=varargin{1};
    end
    filterSpec={'*.xls;*xlsx','Excel files (*.xls,*.xlsx)'};

    [fileName,pathName,~]=uigetfile(...
        filterSpec,'Pick a file','MultiSelect','on');
    
    xlsName=fileName;
    
    oldPath=cd(pathName); % stores old path and changes to new path
    
    [~,sheetNames]=xlsfinfo(fileName);
    nSheets=length(sheetNames);
    
    if nSheets~=1
        [sheetIndex,~]=listdlg('PromptString','Select a file:',...
                    'SelectionMode','single',...
                    'ListString',sheetNames);
        xlsSheet=sheetNames{sheetIndex};
    else
        xlsSheet=sheetNames{1};
    end

    % extract experimental Kp values from Excel file and assign each col to
    % appropriate pressure trials
    [~,~,xlsData] = xlsread(xlsName,xlsSheet);
    iRow=6;
    nRows=size(xlsData,1);
    nCols=size(xlsData,2);

    iCount=1;
    for c=2:nCols
        if ~all(isnan(cell2mat(xlsData(iRow:end,c))))
            dataCols(iCount)=c;
            iCount=iCount+1;
        else
            nCols=nCols-1;
        end
    end
    dataRows=iRow:nRows;
    dataSetName=cell(1,nCols-1);
    
    %parsName={'B1','B2','B3','C1','C2','C3','gamma','lambda'};
    parsName={'B1','B2','B3','C1','C2','C3'};
    nPars=length(parsName);
 
    i=1;
    for c=dataCols
        dataSetName{1,i}=xlsData{2,c};
        i=i+1;
    end
    
    % organize data, variable index 1 of table = col 1 of spreadsheet
    y_exp = cell2mat(xlsData(iRow:end,1)); % store entire yCO2 col
    Kp_exp = zeros(nRows-(iRow-1),nCols-1);
    for iCol=1:nCols-1
        Kp_exp(:,iCol)=cell2mat(xlsData(iRow:end,iCol+1));
    end
    
    % The pure species isotherm values at T,P equaling at infinite dilution
    % of the binaries:
    userVals = cell2mat(xlsData(iRow-2:iRow-1,2:end));
    
    figure; 
    if(nCols-1==1)
        nPlotRows=1;
        nPlotCols=1;
    elseif(nCols-1==2)
        nPlotRows=1;
        nPlotCols=2;
    elseif(nCols-1==3)
        nPlotRows=1;
        nPlotCols=3;    
    elseif(nCols-1==4)
        nPlotRows=2;
        nPlotCols=2;
    elseif(nCols-1==5)
        nPlotRows=2;
        nPlotCols=3;
    elseif(nCols-1>=6)
        nPlotRows=ceil((nCols-1)/3);
        nPlotCols=3;
    end
    
    LB = [-Inf -Inf -Inf -Inf -Inf -Inf];
    UB = Inf.*ones(1,nPars);
    
    
    
    iCount=1;
    for i=1:size(Kp_exp,2)
        
        options=optimoptions(@patternsearch);
        if exist('initialPars','var')
            x0=initialPars(iCount).pars;
        else
            x0=randn(nPars,1);
        end
        
        SSR_fun = residual(y_exp,Kp_exp(:,iCount));
        
        [fitPars_min,SSR_min] = patternsearch(SSR_fun,x0,[],[],[],[],LB,UB,nlcon_fun(y_exp,Kp_exp(:,i),userVals(:,i)),options);
        
        yModel = linspace(0,1,100);
        
        fitPars=fitPars_min;
        SSR=SSR_min;
        
        answer = {'no','yes'};
        
        [~,dqdP,q] = Kp_fun(yModel,fitPars);
        fprintf( ...
  '\nTrial: %s \nSSR=%.4g \n[B1,B2,B3,C1,C2,C3]=[%s] \ndqdP>0: %s \nq(y=0)=%.4g \nq(y=1)=%.4g \n',...
            dataSetName{i},SSR_min,sprintf('%#.4g ',fitPars_min),answer{double(min(dqdP(:))>=0) + 1},q(2,1),q(1,end));

        subplot(nPlotRows,nPlotCols,i);
        plot(y_exp,Kp_exp(:,i),'or',...
            yModel,Kp_fun(yModel,fitPars));
        ylim([0 1.10*max(Kp_exp(1,i))]);
        xlim([0 1]);
        xl=xlim;
        yl=ylim;
        legend('Experimental','Fitted Model','location','southeast');
        title(sprintf('Data Set: %s (SSR=%.4f)',dataSetName{i},SSR_min));
        %text(0.5*xl(2),0.5*yl(2),sprintf('SSR=%.4f',SSR_min));  
        
        % Return ga_fitting return variable as fits (a
        % structure) with fields pressures, y_exp,
        % Kp_exp, SSR, and pars.
        fits(iCount).dataSetName    = dataSetName{i};
        fits(iCount).y_exp          = y_exp;
        fits(iCount).yModel         = yModel;
        fits(iCount).Kp_exp         = Kp_exp;
        fits(iCount).SSR            = SSR;
        fits(iCount).pars           = fitPars;
        fits(iCount).dqdP           = dqdP;
        fits(iCount).q              = q;
        
        iCount=iCount+1;
    end
    
    figure;
    iCount=1;
    for i=1:(nCols-1)
        subplot(nPlotRows,nPlotCols,i);
        
        plot(fits(iCount).yModel,fits(iCount).q);
        
%         ylim([0 1.10*max(fits(iCount).q)]);
        xlim([0 1]);

        legend('Methane','Nitrogen','location','southeast');
        title(sprintf('q: %s (SSR=%.4f)',dataSetName{i},fits(iCount).SSR));
        iCount=iCount+1;
    end
    
    figure;
    iCount=1;
    for i=1:(nCols-1)
        subplot(nPlotRows,nPlotCols,i);
        
        plot(fits(iCount).yModel,fits(iCount).dqdP);
        
%         ylim([0 1.10*max(fits(iCount).dqdP)]);
        xlim([0 1]);

        legend('Methane','Nitrogen','location','southeast');
        title(sprintf('dqdP: %s (SSR=%.4f)',dataSetName{i},fits(iCount).SSR));
        iCount=iCount+1;
    end
    
    % Prepare output cell array xlsOut to export as .csv file
    xlsOut = xlsData;
    
    xlsOut((nRows+1):(nRows+nPars),1)=parsName;   
    iCount=1;
    for c=dataCols
        xlsOut((nRows+1):(nRows+nPars),c)=num2cell(fits(iCount).pars);
        iCount=iCount+1;
    end
    
    xlsOut(end+1,1)={'SSR'};
    iCount=1;
    for c=dataCols
        xlsOut(end,c)=num2cell(fits(iCount).SSR);
        iCount=iCount+1;
    end   
    
    for i=1:size(xlsOut,1)
        for j=1:size(xlsOut,2)
            if(isnan(xlsOut{i,j}))
                xlsOut{i,j}='';
            end
        end
    end
    
    if(input('\nWould you like to export results to a csv file?\nYes (y) or No (n): ','s')=='y')
        xlsOutName=input('\nEnter a current or new .csv filename, e.g. output.csv: ','s');
        csv_export(xlsOutName,xlsOut);
    end
    
    cd(oldPath);
end

% This is the function containing the equation to fit.
% Returns a function handle to be used in ga(...) 
% inside the ga_fitting function.
% This is the function that you can change if you want to use a different
% equation to fit against. Keep the signature the same (i.e. function val =
% Kp_fun(y,args)) but change the assignments of the the parameter to
% args(index_number) accordingly.
% e.g. add "beta=args(9);" after "theta=args(8);" if beta becomes a new 
% equation parameter used in what follows after "val = ".
function [Kp,dqdP,q] = Kp_fun(y,args)
    b0=args(1); b1=args(2); b2=args(3);
    c0=args(4); c1=args(5); c2=args(6);
    
    % dqdP values calculated to ensure positivity
    dqdP1 = b0+b1.*y+b2.*y.^2;
    dqdP2 = c0+c1.*y+c2.*y.^2;
    dqdP  = [dqdP1' dqdP2'];
    
    % q isotherm values calculated to ensure same as pure species at
    % endpoints to represent infinite dilution (userVals array has the
    % experimental pure species values)
    q1 = b0.*y+(1/2).*b1.*y.^2+(1/3).*b2.*y.^3;
    q2 = c0.*(1-y)+(1/2).*c1.*(1-y.^2)+(1/3).*b2.*(1-y.^3);
    q  = [q1(:).'; q2(:).'];
    
    Kp = (1-y).*(dqdP1) + y.*(dqdP2);
end

% This is the other model function using a complicated mix of logarithm
% exponential, and rational expressions:
function [Kp,dqdP,q] = Kp_fun2(y,args)
    b1=args(1);    b2=args(2); b3=args(3);
    c1=args(4);    c2=args(5); c3=args(6);
    gamma=args(7); theta=args(8);

    dqdP1 = b1./abs(y+gamma) + b2.*exp(theta.*y) + b3;
    dqdP2 = c1./abs(y+gamma) + c2.*exp(theta.*y) + c3;
    dqdP  = [dqdP1' dqdP2'];
    
    % q1=Me q2=N2
    q1 = b1.*sign(y+gamma).*log(y+gamma) + (b2./theta).*exp(theta.*y) + b3.*y;
    q2 = c1.*sign(y+gamma).*log(y+gamma) + (c2./theta).*exp(theta.*y) + c3.*y;
    q  = [q1' q2'];
    
    Kp = (1-y).*(dqdP1) + y.*(dqdP2);
    
    if(~isreal(q))
%         disp(args);
        input('\nq NOT REAL :\n');
    end
    if(~isreal(dqdP))
    %         disp(args);
        input('\nq NOT REAL :\n');
    end
    if(~isreal(Kp))
    %         disp(args);
        input('\nq NOT REAL :\n');
    end
    
    
end

% This is the sum of squared residuals function.
% Returns a function handle to be used in ga(...) 
% inside the ga_fitting function.
function SSR_fun = residual(y_exp,Kp_exp)
    Kp_model = @(args)Kp_fun(y_exp,args);
    
    function val = SSR(args)        
        val = sumsqr(Kp_exp-Kp_model(args));
    end
    
    SSR_fun=@SSR;
end

% This is the nonlinear constraint function.
% Returns a function handle to be used in ga(...) 
% inside the ga_fitting function.
function fun = nlcon_fun(y_exp,Kp_exp,q_pure)

    function [c, ceq] = nlcon(args)  
        [Kp_model,dqdP,q] = Kp_fun(y_exp,args);
        
        c   = [ sumsqr(Kp_exp-Kp_model)-0.0000002;
                -abs(min(dqdP(:)))+0.001;
                -abs(min(q(:)))+0.001];  
            
%         ceq = [ Kp_model(1)-Kp_exp(1);
%                 Kp_model(end)-Kp_exp(end) ];
            
        ceq = [ Kp_model(1)-Kp_exp(1);
                Kp_model(end)-Kp_exp(end);
                q(2,1)-q_pure(1);
                q(1,end)-q_pure(2)];
     end
    
    fun=@nlcon;
end

% This is a utility function to allow exporting
% ga_fitting results to a .csv file if the user 
% desires to do so. Taken online from:
% http://www.mathworks.com/matlabcentral/fileexchange/48560-csvexport-filename-cellvals-
function out = csv_export(filename,cellVals)
    % (C) Nov 25, 2014, Pratik Chhatbar, Feng Lab
    % MUSC Stroke Center, Charleston, SC
    % chhatbar@musc.edu, pratikchhatbar@gmail.com
    
    out=0;
    if length(filename)<5 || ~strncmp(filename(end-3:end),'.csv',4)
        outfile=[filename '.csv'];
    else
        outfile=filename;
    end
    
    fh = fopen(outfile,'w'); % open a file with write privileges, will overwrite old versions
    for ii = 1:size(cellVals,1)
        for  jj = 1:size(cellVals,2)
            if jj==1
                addcoma=[];
            else
                addcoma=', ';
            end
            if ischar(cellVals{ii,jj})
                fwrite(fh,[addcoma,cellVals{ii,jj}]);
            else
                fwrite(fh,[addcoma,num2str(cellVals{ii,jj},'%f')]);
            end
        end
        fwrite(fh,sprintf('\r\n')); % print line break
    end
    fclose(fh); % close file out when done writing
    out=1;
end
