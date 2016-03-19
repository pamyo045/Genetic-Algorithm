# genetic-algorithm
## Contents
1. [Functions](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#1-functions)
  1. [ga_fitting.m](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#1i-ga_fittingm)
  2. [csvexport.m](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#1ii-csvexportm)
2. [Usage](https://github.com/pamyo045/genetic-algorithm/blob/master/README.md#2-usage)

***
## 1. Functions
### 1.i. ga_fitting.m
This function uses MATLAB's built-in function 'ga(...)' to perform a genetic algorithm in order to find the minimum of an objective function (in this case the smallest sum of squared residual (SSR) between model parameters and experimental values within an excel file).

**The MATLAB code:**
```
function fits = ga_isotherm
    % XLSname: Binary - N2-CO2 with HISIV3000 silicalite combo graph.xlsx
    % XLSsheet: Kp exp - curve fits
    % XLSrange: A1:F22
    
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
    
    fits.pressures = pressures;
    fits.y_exp = y_exp;
    fits.Kp_exp = Kp_exp;
    fits.SSR = SSR;
    fits.pars = fitPars';

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
        csvexport(xlsOutName,xlsOut);
    end
    
    cd(oldPath);
end

function val = Kp_fun(y1,args)
    b1=args(:,1); b2=args(:,2); b3=args(:,3);
    c1=args(:,4); c2=args(:,5); c3=args(:,6);
    gamma=args(:,7); theta=args(:,8);
    
    val = ((1-y1).*(b1./abs(y1+gamma)+...
        b2.*exp(theta.*y1)+b3)+...
        y1.*(c1./abs(y1+gamma)+c2.*exp(theta.*y1)+c3));
end

function SSR_fun = residual(y_exp,Kp_exp)
    Kp_model = @(args)Kp_fun(y_exp,args);
    
    function val = SSR(x)
        b1=x(:,1); b2=x(:,2); b3=x(:,3);
        c1=x(:,4); c2=x(:,5); c3=x(:,6);
        gamma=x(:,7); theta=x(:,8);
        args=[b1,b2,b3,c1,c2,c3,gamma,theta];
        
        val = sumsqr(Kp_exp-Kp_model(args));
    end
    
    SSR_fun=@SSR;
end

function fun = nlcon_fun(y_exp,Kp_exp)
    Kp_model = @(args)Kp_fun(y_exp,args);
        
    function [c, ceq] = nlcon(x)
        b1=x(:,1); b2=x(:,2); b3=x(:,3);
        c1=x(:,4); c2=x(:,5); c3=x(:,6);
        gamma=x(:,7); theta=x(:,8);
        args=[b1,b2,b3,c1,c2,c3,gamma,theta];
        
        c = [sumsqr(Kp_exp-Kp_model(args))-0.2];
        ceq = [];
    end
    
    fun=@nlcon;
    
end
```

### 1.ii. csvexport.m
This function was taken from http://www.mathworks.com/matlabcentral/fileexchange/48560-csvexport-filename-cellvals-/content//csvexport.m in order to be able to export results to a .csv file if the user desires to do so and is compatible for all OS platforms (MATLAB's built-in csvwrite and xlswrite have difficulties running on MacOS platforms).
***
## 2. Usage
* Run ga_isotherm.m by either having ga_isotherm.m open in the MATLAB editor and pressing on the 'Run' button, or
  Enter "ga_isotherm" in MATLAB's command window and press Enter.
* Follow the instructions to select the appropriate file path, sheet, and cell range of the Excel file containing the experimental data
  * Note: the range of the data must include the layout shown in **Fig. 1**. Make sure to respect this template style but feel free to have any number of columns (only one column of Kp values are required). For example, for **Fig. 1** the proper input for the range when prompt would be to type "A1:F22" without the quotes.
![fig1](https://github.com/pamyo045/genetic-algorithm/blob/master/Resources/Excel%20Input%20Data%20Template.png)
**Fig. 1:** Excel Input Data Template
* When prompt, decide on exporting the results to a .csv file (either already existing in the path folder of the ga_isotherm.m or enter a new name to create a new file).
  * e.g. type "results.csv" without the quotes and press Enter. This will either overwrite a file named result.csv if it already exists or create a new one if not.
