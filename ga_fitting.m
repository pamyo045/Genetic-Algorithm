function fits = ga_fitting
    % xlsName: Binary - N2-CO2 with HISIV3000 silicalite combo graph.xlsx
    % xlsSheet: Kp exp - curve fits
    % xlsRange: A1:F22
    
    filterSpec={'*.xls;*xlsx','Excel files (*.xls,*.xlsx)'};

    [fileName,pathName,~]=uigetfile(...
        filterSpec,'Pick a file','MultiSelect','on');
    
    xlsName=fileName;
    
    oldPath=cd(pathName); % stores old path and changes to new path
    
    [~,sheetNames]=xlsfinfo(fileName);
    nSheets=length(sheetNames);
    
    [sheetIndex,~]=listdlg('PromptString','Select a file:',...
                'SelectionMode','single',...
                'ListString',sheetNames);
    xlsSheet=sheetNames{sheetIndex};
    
    prompt = 'Enter Excel sheet range of the experimental data:';
    dlgTitle = 'Select range';
    numLines = 1;
    defRange = {'A1:F22'};
    xlsRange = inputdlg(prompt,dlgTitle,numLines,defRange);

    % extract experimental Kp values from Excel file and assign each col to
    % appropriate pressure trials
    [~,~,xlsData] = xlsread(xlsName,xlsSheet,xlsRange{1});
    
    nRows=size(xlsData,1);
    nCols=size(xlsData,2);
    pressures=zeros(1,nCols-1);
    
    nPars=8;
    parsName={'B1','B2','B3','C1','C2','C3','gamma','lambda'};
    for c=1:nCols-1
        pressures(1,c)=str2double(cell2mat(regexp(xlsData{2,1+c},'\d*','match')));
    end
    
    iRow=4;
    % organize data, variable index 1 of table = col 1 of spreadsheet
    y_exp = cell2mat(xlsData(iRow:end,1)); % store entire yCO2 col
    Kp_exp = zeros(nRows-(iRow-1),nCols-1);
    for iCol=1:nCols-1
        Kp_exp(:,iCol)=cell2mat(xlsData(iRow:end,iCol+1));
    end

    figure;
    fitPars=zeros(nCols-1,nPars);
    SSR=zeros(1,nCols-1);
    for i=1:nCols-1
        SSR_fun = residual(y_exp,Kp_exp(:,i));

        options = gaoptimset('vectorized','off','TolFun',1e-9);
        [fitPars_min,SSR_min] = ga(SSR_fun,nPars,[],[],[],[],[],[],nlcon_fun(y_exp,Kp_exp(:,i)),options);
        fprintf('Trial: %d atm \nSSR=%.4g \n[B1,B2,B3,C1,C2,C3,gamma,lambda]=[%s]\n',...
            pressures(i),SSR_min,sprintf('%#.4g ',fitPars_min));
              
        fitPars(i,:)=fitPars_min(:);
        SSR(1,i)=SSR_min;
        yModel = linspace(0,1,100);
        
        subplot(3,2,i);
        plot(y_exp,Kp_exp(:,i),'or',...
            yModel,Kp_fun(yModel,fitPars(i,:))),'-b';
        ylim([0 1.10*Kp_exp(1,i)]);
        xlim([0 1]);
        xl=xlim;
        yl=ylim;
        legend('Experimental','Fitted Model');
        title(sprintf('System Pressure: %d atm',i));
        text(0.5*xl(2),0.5*yl(2),sprintf('SSR=%.4f',SSR_min));            
    end
    
    % Return ga_fitting return variable as fits (a 
    % structure) with fields pressures, y_exp,
    % Kp_exp, SSR, and pars.
    fits.pressures = pressures;
    fits.y_exp = y_exp;
    fits.Kp_exp = Kp_exp;
    fits.SSR = SSR;
    fits.pars = fitPars';

    % Prepare output cell array xlsOut to export as .csv file
    xlsOut = xlsData;
    xlsOut((nRows+1):(nRows+nPars),1)=parsName;
    xlsOut((nRows+1):(nRows+nPars),2:nCols)=num2cell(fits.pars);
    xlsOut(end+1,1)={'SSR'};
    xlsOut(end,2:nCols)=num2cell(fits.SSR);
    for i=1:size(xlsOut,1)
        for j=1:size(xlsOut,2)
            if(isnan(xlsOut{i,j}))
                xlsOut{i,j}='';
            end
        end
    end
    
    if(input('Would you like to export results to a csv file?\nYes (y) or No (n): ','s')=='y')
        xlsOutName=input('Enter a current or new .csv filename, e.g. output.csv: ','s');
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
function val = Kp_fun(y,args)
    b1=args(1);    b2=args(2); b3=args(3);
    c1=args(4);    c2=args(5); c3=args(6);
    gamma=args(7); theta=args(8);
    
    val = ((1-y).*(b1./abs(y+gamma)+...
        b2.*exp(theta.*y)+b3)+...
        y.*(c1./abs(y+gamma)+c2.*exp(theta.*y)+c3));
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
function fun = nlcon_fun(y_exp,Kp_exp)
    Kp_model = @(args)Kp_fun(y_exp,args);
        
    function [c, ceq] = nlcon(args)     
        c = [sumsqr(Kp_exp-Kp_model(args))-0.2];
        ceq = [];
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
