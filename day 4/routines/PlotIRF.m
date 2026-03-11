function PlotIRF(IRFinput,COLORPICK,LABELTITLE,YlabelTitle)

if nargin == 3 
    YlabelTitle = ['\%  Dev.'];
end

Shade(1) = .92; 
Shade(2) = .86; 
switch COLORPICK
    case 'grey' 
        ColorLine = 'b'; % If Grey Use Blue color
        ShadeIntensity = .85*[1 1 1]; 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2); 
        LineStyleSet = '-';
    case 'blue' 
        ColorLine = 'b';
        ShadeIntensity = [0.5 .8 1]; 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2); 
        LineStyleSet = '-';
    case 'red'
        ColorLine = 'r'; 
        ShadeIntensity = ([252 211 210]/256); 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2);  
        LineStyleSet = ':';
    case 'green'
        ColorLine = [0.45 0.9 0.15]; %'g'; 
        ShadeIntensity = ([0.84 0.91 0.85]); 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2);  
        LineStyleSet = '-';
    case 'magenta'
        ColorLine = 'm'; %[1 0 1]; 
        ShadeIntensity = ([1 0.87 1]); 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2);  
        LineStyleSet = ':';    
    case 'black' 
        ColorLine = 'k'; % If Grey Use Blue color
        ShadeIntensity = .85*[1 1 1]; 
        ShadeIntensity68 = ShadeIntensity;
        ShadeIntensity = ShadeIntensity*Shade(1); 
        ShadeIntensity68 = ShadeIntensity68*Shade(2); 
        LineStyleSet = '-';    
end

        


%--------------------------------------------------------------------------
%             PLOT FEATURES
%--------------------------------------------------------------------------

FontS = 12+5; % Setting the fontsize for the axes scales
FontT = 14+6; % Setting the fontsize for the title
fontAx = 12+5; % setting fontsize of axes labels
LineType =  LineStyleSet; %'-ob';     % line formatting
LineW = 2+2; % setting the linewidth
MarkSize = 5; % marker size
INTERPRETERPICK = 'latex'; 

% extEPS = '.eps'; % for .eps format
% extPDF = '.pdf'; % for .pdf format
% extFIG = '.fig'; % for .fig format
% % Save = 1;       % set to 1 if you want MATLAB to save automatically after plotting
% % format = 1;     % =1 if .eps; =2 if pdf; =3 if .fig;

%--------------------------------------------------------------------------
%             PRE-PROCESSING DATA 
%--------------------------------------------------------------------------
% 1. IRFs
irfhorizon = size(IRFinput,1); 
figirfhorizon = irfhorizon;
MeanIRF = IRFinput(:,3)'; 
% 2. Confidence Bands
IRFLB90 = IRFinput(:,1)';
IRFLB68 = IRFinput(:,2)';
IRFUB68 = IRFinput(:,4)';
IRFUB90 = IRFinput(:,5)';

% 3. Zero line
zeroline = zeros(irfhorizon,1);
x = [0:irfhorizon-1]';

%--------------------------------------------------------------------------
%            PLOT SERIES and ERROR BANDS
%--------------------------------------------------------------------------
% 1. Confidence bands / shade
ha = shadedplot(x, IRFUB90,IRFLB90, ShadeIntensity, ShadeIntensity);
hold on;
hb = shadedplot(x, IRFUB68,IRFLB68, ShadeIntensity68, ShadeIntensity68);
hold on;
% 2. IRF
plot(x, MeanIRF,'linestyle',LineType,'color',ColorLine,'MarkerSize',MarkSize,'LineWidth',LineW);
hold on;
% Zero line
plot(zeroline(1:figirfhorizon),'k')
box on
grid on;
hold off
%--------------------------------------------------------------------------
%                               PLOT FORMATTING
%--------------------------------------------------------------------------
% Y axis Formatting
ylabel(YlabelTitle,'fontsize',fontAx,'Interpreter',INTERPRETERPICK);

% X axis formatting
set(gca,'XTick',[0:1:figirfhorizon-1],'XLim',[0 figirfhorizon-1],'fontsize',FontS) % set x-axis span
xlabel('Years','fontsize',fontAx,'Interpreter',INTERPRETERPICK);

% Titles
title(LABELTITLE,'fontweight','bold','fontsize',FontT,'Interpreter',INTERPRETERPICK)
% title(LABELTITLE,'fontsize',FontT,'Interpreter',INTERPRETERPICK)
        
%--------------------------------------------------------------------------
% Plotting Area
box on;
set(gcf, 'Color', 'w');
%--------------------------------------------------------------------------

