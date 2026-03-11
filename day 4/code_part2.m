%% Question 4: Panel LP
clear all
close all 
clc

addpath("routines");
addpath(genpath("PanelDataToolboxMatlab"));
addpath("mat_files");
addpath("panel_LP");

%% LOAD DATA 
load("SaveData_IDB_Global.mat"); % Dollar index
load("SaveData_IDB_Final.mat"); % Dates, country labels, data
load("Store_InstrumentRisk_Short.mat"); % Instruments for uncertainty

%% SETUP LP 
settingsLP.maxH = 5; % Maximum number of horizons
settingsLP.MaxLPLagsOwn = 2; % Maximum number of lags of dependent variable
settingsLP.MaxLPLagsOther = 2; % Maximum number of lags of other controls

%% DECLARE CASES GDP, PX, TB, inflation, PB, RER
WhichVariables{1,1} = 'GDP'; 
WhichVariables{size(WhichVariables,1)+1,1} = 'Export Prices' ;
WhichVariables{size(WhichVariables,1)+1,1} = 'Trade Balance' ;
WhichVariables{size(WhichVariables,1)+1,1} = 'Inflation' ;
WhichVariables{size(WhichVariables,1)+1,1} = 'Primary Balance' ;
WhichVariables{size(WhichVariables,1)+1,1} = 'RER' ;


%% Orgnize outputs + define common inputs
No_Variables = size(WhichVariables,1); % Number of LPs, i.e. number of dependent variables
No_Countries = size(CountryNames,1); % Number of countries

Collect.Coeff = ones(settingsLP.maxH+1,No_Variables); % (horizon+1) x (number of LPs)
Collect.StdErr = zeros(settingsLP.maxH+1,No_Variables); % Same

Tt = size(FinalData.Px_Detr,1); % Time sample length
Data.id = repmat(1:No_Countries,Tt,1);  clear Tt % (time sample) x (number of countries)
% --> unit identifier
Data.year = repmat(Dates,1,No_Countries); % Same
% --> year identifier

%% Compute panel-LP-IV
for PickVariable = 1:No_Variables % Loop over variables

    XX = DataGlobal.LogDollarIndex; % Variable to instrument --> ΔXX inside the function

    InstrumentAll = [InstrumentRisk.InstrumentGOLD InstrumentRisk.InstrumentFinUnc];  
    % Set of instruments 

    % cat(3,...) concatenates array along their third dimension
    Controls = cat(3,MatrixDiff(FinalData.Px_Detr)); % Log Px                            
    Controls = cat(3,Controls,MatrixDiff(FinalData.GDP)); % Log GDP
    Controls = cat(3,Controls,MatrixDiff(repmat(DataGlobal.LogDollarIndex,1,No_Countries)));
    % repmat(...) replicates the dollar index for the every country in the
    % list, because it has common values across countries
    Controls = cat(3,Controls,MatrixDiff(FinalData.RER)); % RER

    ControlsNames{1} = 'Px change'; 
    ControlsNames{2} = 'GDP growth'; 
    ControlsNames{3} = 'Dollar index'; 
    ControlsNames{4} = 'RER'; 
    
    chooseVariable = WhichVariables{PickVariable,1}; 

    switch chooseVariable 
        case 'GDP'
            yy = FinalData.GDP; 
            ControlsNames = {ControlsNames{[1,3:size(Controls,3)]}}; 
            Controls = Controls(:,:,[1,3:size(Controls,3)]); % Select controls      
            YlabelTitle = '\%  Dev. from s.s.'; 

        case 'Export Prices'  
            yy = FinalData.Px_Detr; 
            Controls = Controls(:,:,3);
            ControlsNames = {ControlsNames{3}};               
            YlabelTitle = '\%  Dev. from s.s.';

        case 'Trade Balance'
            yy = FinalData.TB; 
            ControlsNames = {ControlsNames{[1,2:size(Controls,3)]}}; 
            Controls = Controls(:,:,[1,2:size(Controls,3)]);        
            YlabelTitle = '\%  of GDP'; 

        case 'RER'
            yy = FinalData.RER; 
            ControlsNames = {ControlsNames{[1,2:size(Controls,3)-1]}}; 
            Controls = Controls(:,:,[1,2:size(Controls,3)-1]);                
            YlabelTitle = '\%  Dev. from s.s.'; 
  
        case 'Primary Balance'
            yy = (FinalData.PB); 
            ControlsNames = {ControlsNames{[1,2:size(Controls,3)]}}; 
            Controls = Controls(:,:,[1,2:size(Controls,3)]);        
            YlabelTitle = '\%  of GDP'; 

        case 'Inflation'
            yy = FinalData.LogPrice; 
            ControlsNames = {ControlsNames{[1,2:size(Controls,3)]}}; 
            Controls = Controls(:,:,[1,2:size(Controls,3)]);       
            YlabelTitle = '\%';      
	    end
    
    
    [Collect_LP]=Panel_LP_IV(yy,XX,Controls,[],InstrumentAll,Data.id,Data.year,settingsLP);
    
    Collect.Coeff(:,PickVariable) = Collect_LP.Coeff; 
    Collect.StdErr(:,PickVariable) = Collect_LP.StdErr; 

%% PLOT

        BandsPlot1 = Collect.Coeff(:,PickVariable)*ones(1,2) + Collect.StdErr(:,PickVariable)*[norminv((1-.68)/2) -norminv((1-.68)/2)];
        BandsPlot2 = Collect.Coeff(:,PickVariable)*ones(1,2) + Collect.StdErr(:,PickVariable)*[norminv((1-.9)/2) -norminv((1-.9)/2)];
        MeanPlot = Collect.Coeff(:,PickVariable); 
        IRFinput(:,[1,5,2,4,3]) = [BandsPlot2 BandsPlot1 MeanPlot];
        PlotTitle = [WhichVariables{PickVariable,1}]; 
        PlotColor = 'blue';    

        h = figure; 
        PlotIRF(IRFinput,PlotColor,PlotTitle,YlabelTitle)

        % Ensure zero is included in y-axis
        y_limits = ylim;
        ylim([min(y_limits(1), 0), max(y_limits(2), 0)])
        
                 
end
