function MinimalGlacierModel

% Minimal glacier model Oerlemans

close all;
clear all;

f = figure('Name','Minimal Glacier Model','numbertitle','off', ...
    'Position', [200, 100, 900, 600]);

% Constants
Er0=-200;
slope0=0.06;
gamma=-0.007;
tau0=80e3;
rho=900;
g=9.81;
am=3;
mu=10;
Hm0=100;

% Initialization
Er=Er0;
slope=slope0;
L0=10e3; % initial glacier length
hysteresis=0;
linear=1;


% Other control
npts=500;
alpha=linspace(0.03,0.15,npts);
Ls=NaN(npts,1);
Ls2=NaN(npts,1);
x=linspace(0,30e3,npts)';
Ers=linspace(-400,200,npts)';
ELA=x*0+Er;
b=-slope*x;

%// Make initial plot
[L,h,Hm] = GlacierModel(Hm0,linear,slope,Er,am,mu,tau0,rho,g,x);
p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);

%// re-position the axes to make room for the sliders/buttons
set(gca, 'position', [0.1 0.25 0.8 0.7]);


%// initialize the sliders/buttons
h1 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.05 0.05 0.2 0.05],...
    'min'     , -400,...               %// Make the A between 1...
    'max'     , 200,...              %// and 10, with initial value
    'value'   , Er,...               %// as set above.
    'callback', @sliderEr);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h2 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.28 0.05 0.2 0.05],...
    'min'     , 0.03,...               %// Make the A between 1...
    'max'     , 0.15,...              %// and 10, with initial value
    'value'   , slope,...               %// as set above.
    'callback', @sliderSlope);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
                                    
h3 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...        
    'position', [0.51 0.05 0.1 0.05],...
    'String'  , 'Linear|Nonlinear',...
    'callback', @popupLin);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h4 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...    
    'String'  , 'Glacier|Sensitivity slope|Sensitivity Er',...
    'position', [0.65 0.05 0.10 0.05],...
    'callback', @popupHyst);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar

h = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'pushbutton',...    
    'String'  , 'INIT',...
    'position', [0.8 0.05 0.08 0.05],...
    'callback', @pushInit);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
%// THE MAGIC INGREDIENT
%// ===========================

% hLstn = handle.listener(h1,'ActionEvent',@sliderX0); %#ok<NASGU>
% hLstn = handle.listener(h2,'ActionEvent',@sliderBeta); %#ok<NASGU>
% hLstn = handle.listener(h3,'ActionEvent',@popupLin); %#ok<NASGU>
% hLstn = handle.listener(h4,'ActionEvent',@popupHyst); %#ok<NASGU>
% hLstn = handle.listener(h,'ActionEvent',@pushInit); %#ok<NASGU>

%// (variable appears unused, but not assigning it to anything means that 
%// the listener is stored in the 'ans' variable. If "ans" is overwritten, 
%// the listener goes out of scope and is thus destroyed, and thus, it no 
%// longer works.

%// ===========================


%// The slider's callback:
%//    1) clears the old plot
%//    2) computes new values using the (continuously) updated slider values
%//    3) re-draw the plot and re-set the axes settings

function sliderEr(~,~)
    Er=get(h1,'value');
    ELA=x*0+Er;
    [L,h,Hm] = GlacierModel(Hm,linear,slope,Er,am,mu,tau0,rho,g,x);
    if hysteresis==0
        p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);
    else
        if hysteresis==1
            [Ls]=SensitivitySlope(Hm0,linear,alpha,Er,am,mu,tau0,rho,g,npts);
            p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma);
        else
            [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Ers,am,mu,tau0,rho,g,npts);
            p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma);
        end
    end
end

function sliderSlope(~,~)
    slope=get(h2,'value');
    b=-slope*x;
    [L,h,Hm] = GlacierModel(Hm,linear,slope,Er,am,mu,tau0,rho,g,x);
    if hysteresis==0
        p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);
    else
        if hysteresis==1
            [Ls]=SensitivitySlope(Hm0,linear,alpha,Er,am,mu,tau0,rho,g,npts);
            p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma);
        else
            [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Ers,am,mu,tau0,rho,g,npts);
            p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma);
        end
    end
end

function popupLin(~,~)
    getLin=get(h3,'value');
    if getLin==1
        linear=1;
    else
        linear=0;
    end
    [L,h,Hm] = GlacierModel(Hm,linear,slope,Er,am,mu,tau0,rho,g,x);
    if hysteresis==0
        p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);
    else
        if hysteresis==1
            [Ls]=SensitivitySlope(Hm0,linear,alpha,Er,am,mu,tau0,rho,g,npts);
            p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma);
        else
            [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Ers,am,mu,tau0,rho,g,npts);
            p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma);
        end
    end
end

function pushInit(~,~)
    delete(p);
    set(h1,'Value',Er0);
    set(h2,'Value',slope0);
    Er=Er0;
    slope=slope0;
    L=L0;
    ELA=x*0+Er;
    b=-slope*x;
    [L,h,Hm] = GlacierModel(Hm,linear,slope,Er,am,mu,tau0,rho,g,x);
    if hysteresis==0
        p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);
    else
        if hysteresis==1
            [Ls]=SensitivitySlope(Hm0,linear,alpha,Er,am,mu,tau0,rho,g,npts);
            p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma);
        else
            [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Ers,am,mu,tau0,rho,g,npts);
            p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma);
        end
    end
end

function popupHyst(~,~)
    delete(p);
    graphtype=get(h4,'value');
    if graphtype==1
        hysteresis=0;
        [L,h,Hm] = GlacierModel(Hm,linear,slope,Er,am,mu,tau0,rho,g,x);
        p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma);
    else
        if graphtype==2
            hysteresis=1;
            [Ls]=SensitivitySlope(Hm0,linear,alpha,Er,am,mu,tau0,rho,g,npts);
            p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma);
        else
            hysteresis=2;
            [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Ers,am,mu,tau0,rho,g,npts);
            p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma);
        end
    end
end

end


function [L,h,Hm] = GlacierModel(Hm0,linear,slope,Er,am,mu,tau0,rho,g,x)

if linear==1
% linear model
    Hm=tau0/(rho*g*slope);
    L=2*(Hm-Er)/slope;
    L(L<0)=0;
else
% nonlinear model
    Dis=(2*am/(slope*(1+mu*slope)))^2-8*Er/slope;
    L=0;
    L(Dis>0)=(am/(slope*(1+mu*slope))+0.5*sqrt(Dis))^2;
    L(Er-Hm0>0)=0;
    Hm=am*sqrt(L)/(1+mu*slope);
end

if L>0
    sigma=(1.5*Hm)^2/L;
    h=sqrt(sigma*(L-x));
    h(L-x<0)=0;
else
    h=x*0;
end

end


function p=PlotGlacier(Er,Hm,L,h,x,b,ELA,slope,Er0,gamma)

deltaT=(Er0-Er)*gamma;

p=plot(x*1e-3,b+h,'Linewidth',3); hold on;
plot(x*1e-3,b,'k','Linewidth',3);
plot(x*1e-3,ELA,'r--','Linewidth',3);
if L>=0
    plot(L*1e-3,-slope*L,'o','Color','k','MarkerFaceColor','k', ...
        'MarkerSize',12);
end
grid on;
axis([0 25 -1500 500]);
xlabel('Horizontal distance (km)');
ylabel('Elevation (m from b_0)');

plotval=strcat('E_r = ',num2str(Er));
text(0,-1850,plotval);
plotval=strcat('s = ',num2str(slope));
text(8,-1850,plotval);
plotval=strcat('L = ',num2str(L*1e-3),' km');
text(17,-1750,plotval);
plotval=strcat('H_m = ',num2str(Hm),' m');
text(22,-1750,plotval);
plotval=strcat('\Delta T = ',num2str(deltaT),' °C');
text(1,-1400,plotval);
text(20,400,'\bf \fontsize{12} Glacier surface','BackgroundColor','w','Color','b');
text(20,300,'\bf \fontsize{12} Bedrock','BackgroundColor','w','Color','k');
text(20,200,'\bf \fontsize{12} Equilibrium line','BackgroundColor','w','Color','r');
hold off;

end


function [Ls]=SensitivitySlope(Hm0,linear,slope,Er,am,mu,tau0,rho,g,npts)

if linear==1
% linear model
    Hm=tau0./(rho*g*slope);
    Ls=2*(Hm-Er)./slope;
else
% nonlinear model
    Ls=NaN(npts,1);
    Dis=(2*am./(slope.*(1+mu*slope))).^2-8*Er./slope;
    Ls(Dis>0)=(am./(slope(Dis>0).*(1+mu*slope(Dis>0)))+0.5*sqrt(Dis(Dis>0))).^2;
    Ls(Er-Hm0>0)=0;
end

end


function p=PlotSlope(Ls,alpha,L,slope,Er,Er0,Hm,gamma)

deltaT=(Er0-Er)*gamma;

p=plot(alpha,Ls*1e-3,'Linewidth',3); hold on;
if L>=0
    plot(slope,L*1e-3,'o','Color','k','MarkerFaceColor','k', ...
        'MarkerSize',12);
end
xlabel('Bedrock slope (rad)');
ylabel('Glacier length (km)');
axis([0.02 0.16 0 35]);
grid on;

plotval=strcat('E_r = ',num2str(Er));
text(0.02,-6,plotval);
plotval=strcat('s = ',num2str(slope));
text(0.065,-6,plotval);
plotval=strcat('L = ',num2str(L*1e-3),' km');
text(0.11,-4,plotval);
plotval=strcat('H_m = ',num2str(Hm),' m');
text(0.14,-4,plotval);
plotval=strcat('\Delta T = ',num2str(deltaT),' °C');
text(0.022,2,plotval);

hold off;

end


function [Ls,Ls2]=SensitivityEr(Hm0,linear,slope,Er,am,mu,tau0,rho,g,npts)

if linear==1
% linear model
    Hm=tau0/(rho*g*slope);
    Ls=2*(Hm-Er)/slope;
    Ls2=Ls;
else
% nonlinear model
    Ls=NaN(npts,1);
    Ls2=NaN(npts,1);
    Dis=(2*am/(slope*(1+mu*slope)))^2-8*Er/slope;
    Ls(Dis>0)=(am/(slope*(1+mu*slope))+0.5*sqrt(Dis(Dis>0))).^2;
    Ls2(Dis>0)=(am/(slope*(1+mu*slope))-0.5*sqrt(Dis(Dis>0))).^2;
end

end


function p=PlotEr(Ls,Ls2,Ers,L,slope,Er,Er0,Hm,gamma)

deltaT=(Er0-Er)*gamma;

p=plot(Ers,Ls2*1e-3,'r','Linewidth',3); hold on;
plot(Ers,Ls*1e-3,'Linewidth',3);
if L>=0
    plot(Er,L*1e-3,'o','Color','k','MarkerFaceColor','k', ...
        'MarkerSize',12);
end
xlabel('E_r = E-b_0 (m)');
ylabel('Glacier length (km)');
axis([-400 200 0 35]);
grid on;

plotval=strcat('E_r = ',num2str(Er));
text(-400,-6,plotval);
plotval=strcat('s = ',num2str(slope));
text(-200,-6,plotval);
plotval=strcat('L = ',num2str(L*1e-3),' km');
text(0,-4,plotval);
plotval=strcat('H_m = ',num2str(Hm),' m');
text(100,-4,plotval);
plotval=strcat('\Delta T = ',num2str(deltaT),' °C');
text(-380,2,plotval);
text(80,33,'\bf \fontsize{12} Stable branch','BackgroundColor','w','Color','b');
text(80,31,'\bf \fontsize{12} Unstable branch','BackgroundColor','w','Color','r');

hold off;

end



