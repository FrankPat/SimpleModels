function WeertmanIceSheetModel

% Weertman ice sheet model (1974)

close all;
clear all;

f = figure('Name','Weertman Ice Sheet Model','numbertitle','off', ...
    'Position', [200, 100, 900, 600]);

% Constants
lambda=14;
rho=900;
rhom=3000;
zeta=rho/(rhom-rho);
g=9.81;
% tau0=0.5*lambda*rho*g*(1+zeta);

% Initialization
L0=400;
x00=-400;
beta0=0.002;
Vac0=1.2;
epsilon0=0.24;
hysteresis=0;

% Parameters
L=L0;
x0=x00;
beta=beta0;
Vac=Vac0;
epsilon=epsilon0;
xf=0;
hf=0;


% Other control
delta=20e3;
x=(-5000e3:delta:5000e3)';
xw=linspace(-800e3,200e3,1000);

% Controls
Xmin=-5000;
Xmax=5000;
Ymin=-4000;
Ymax=8000;

%// Make initial plot
[L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x00*1e3,lambda,L0*1e3,x);
p=plotWeertman(L,x,h,hs,xf,hf,x00*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);

%// re-position the axes to make room for the sliders/buttons
set(gca, 'position', [0.1 0.35 0.8 0.6]);

%// initialize the sliders/buttons
h1 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.05 0.05 0.2 0.05],...
    'min'     , -800,...               %// Make the A between 1...
    'max'     , 200,...              %// and 10, with initial value
    'value'   , x0,...               %// as set above.
    'callback', @sliderX0);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h2 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.28 0.05 0.2 0.05],...
    'min'     , 0.0015,...               %// Make the A between 1...
    'max'     , 0.003,...              %// and 10, with initial value
    'value'   , beta,...               %// as set above.
    'callback', @sliderBeta);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
                                    
h3 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.51 0.05 0.2 0.05],...
    'min'     , 0.1,...               %// Make the A between 1...
    'max'     , 0.3,...              %// and 10, with initial value
    'value'   , epsilon,...               %// as set above.
    'callback', @sliderEps);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h4 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...    
    'String'  , 'Ice sheet|Equilibria',...
    'position', [0.8 0.15 0.10 0.05],...
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
% hLstn = handle.listener(h3,'ActionEvent',@sliderEps); %#ok<NASGU>
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

function sliderX0(~,~)
    x0=get(h1,'value');
    [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0*1e3,lambda,L*1e3,x);
    if hysteresis==0
        p=plotWeertman(L,x,h,hs,xf,hf,x0*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);
    else
        [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,xw,lambda);
        p=plotHysteresis(L,hf,xw*1e-3,x0,beta,epsilon,L1,L2, ...
            min(xw*1e-3),max(xw*1e-3),0,Xmax);
    end
end

function sliderBeta(~,~)
    beta=get(h2,'value');
    [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0*1e3,lambda,L*1e3,x);
    if hysteresis==0
        p=plotWeertman(L,x,h,hs,xf,hf,x0*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);
    else
        [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,xw,lambda);
        p=plotHysteresis(L,hf,xw*1e-3,x0,beta,epsilon,L1,L2, ...
            min(xw*1e-3),max(xw*1e-3),0,Xmax);
    end
end

function sliderEps(~,~)
    epsilon=get(h3,'value');
    [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0*1e3,lambda,L*1e3,x);
    if hysteresis==0
        p=plotWeertman(L,x,h,hs,xf,hf,x0*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);
    else
        [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,xw,lambda);
        p=plotHysteresis(L,hf,xw*1e-3,x0,beta,epsilon,L1,L2, ...
            min(xw*1e-3),max(xw*1e-3),0,Xmax);
    end
end

function pushInit(~,~)
    delete(p);
    set(h1,'Value',x00);
    set(h2,'Value',beta0);
    set(h3,'Value',epsilon0);
    x0=x00;
    beta=beta0;
    epsilon=epsilon0;
    L=L0;
    [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0*1e3,lambda,L*1e3,x);
    if hysteresis==0
        p=plotWeertman(L,x,h,hs,xf,hf,x0*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);
    else
        [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,xw,lambda);
        p=plotHysteresis(L,hf,xw*1e-3,x0,beta,epsilon,L1,L2, ...
            min(xw*1e-3),max(xw*1e-3),0,Xmax);
    end
end

function popupHyst(~,~)
    delete(p);
    graphtype=get(h4,'value');
    if graphtype==1
        hysteresis=0;
    else
        hysteresis=1;
    end
    [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0*1e3,lambda,L*1e3,x);
    if hysteresis==0
        p=plotWeertman(L,x,h,hs,xf,hf,x0*1e3,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta);
    else
        [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,xw,lambda);
        p=plotHysteresis(L,hf,xw*1e-3,x0,beta,epsilon,L1,L2, ...
            min(xw*1e-3),max(xw*1e-3),0,Xmax);
    end
end

end



function [L,h,hs,xf,hf]=WeertmanModel(Vac,epsilon,beta,x0,lambda,L,x)

% Initial determination of viability of ice sheet
% based on previous determination of L
h=sqrt(lambda*(L-abs(x)));
h(L-abs(x)<0)=0;

% ELA and position on initial ice sheet
a=beta^2;
b=lambda-2*beta^2*x0;
c=beta^2*x0^2-lambda*L;
Dis=b^2-4*a*c;
xf=(-b+sqrt(Dis))/(2*a);
hs=beta*(x-x0);

h0=sqrt(lambda*L);
hs0=-beta*x0;

if h0>hs0
    % Analytical solution of ice sheet profile
    Vab=Vac/epsilon;
    mu=(Vac+Vab)/(Vab-Vac);
    a=beta^2;
    b=lambda*(1-mu)-2*x0*beta^2;
    c=beta^2*x0^2;
    Dis=b^2-4*a*c;
    if Dis<0
        L=0;
    else
        L=mu*(-b+sqrt(Dis))/(2*a);
    end
    h=sqrt(lambda*(L-abs(x)));
    h(L-abs(x)<0)=0;

    % ELA and position on ice sheet
    a=beta^2;
    b=lambda-2*beta^2*x0;
    c=beta^2*x0^2-lambda*L;
    Dis=b^2-4*a*c;
    xf=(-b+sqrt(Dis))/(2*a);
    if xf>=0
        hf=sqrt(lambda*(L-abs(xf)));
    else
        hf=0;
        xf=0;
    end
else
    L=0;
    xf=0;
    hf=0;
end
L=L*1e-3;

end

function p=plotWeertman(L,x,h,hs,xf,hf,x0,beta,epsilon,Xmin,Xmax,Ymin,Ymax,zeta)

p=plot(x*1e-3,h,'Linewidth',3); hold on;
plot(x*1e-3,hs,'r--','Linewidth',3);
plot(x*1e-3,-zeta*h,'k','Linewidth',3);
plot(xf*1e-3,hf,'ok','Linewidth',4,'Color','k','MarkerFaceColor','k', ...
    'MarkerSize',12);
plot(x0*1e-3,0,'or','Linewidth',3,'Color','r');
axis([Xmin Xmax Ymin Ymax]);
grid on;
plotval=strcat('x_0 = ',num2str(x0*1e-3));
text(Xmin,Ymin-4500,plotval);
plotval=strcat('\beta = ',num2str(beta));
text(-2000,Ymin-4500,plotval);
plotval=strcat('\epsilon = ',num2str(epsilon));
text(1000,Ymin-4500,plotval);
plotval=strcat('L = ',num2str(L),' km');
text(-2000,Ymin-3000,plotval);
plotval=strcat('ELA = ',num2str(hf),' m a.s.l.');
text(0,Ymin-3000,plotval);
hold off;
xlabel('Horizontal distance from center (km)');
ylabel('Elevation (m a.s.l.)');
text(-4900,7500,'\bf \fontsize{12} Ice sheet surface','BackgroundColor','w','Color','b');
text(-4900,6700,'\bf \fontsize{12} Bedrock','BackgroundColor','w','Color','k');
text(-4900,5900,'\bf \fontsize{12} Equilibrium line','BackgroundColor','w','Color','r');

end

function [L1,L2]=WeertmanLfunc(Vac,epsilon,beta,x0,lambda)

% Makes hysteresis plot

Vab=Vac/epsilon;
L1=NaN(length(x0),1);
L2=L1;
% Analytical solution
mu=(Vac+Vab)/(Vab-Vac);
a=beta^2;
b=lambda*(1-mu)-2*x0*beta^2;
c=beta^2*x0.^2;
Dis=b.^2-4*a*c;
L1(Dis>=0)=mu*(-b(Dis>=0)+sqrt(Dis(Dis>=0)))/(2*a);
L2(Dis>=0)=mu*(-b(Dis>=0)-sqrt(Dis(Dis>=0)))/(2*a);
L1=L1*1e-3;
L2=L2*1e-3;

end

function p=plotHysteresis(L,hf,xw,x0,beta,epsilon,L1,L2,Xmin,Xmax,Ymin,Ymax)

p=plot(xw,L1,'b','Linewidth',3); hold on;
plot(xw,L2,'r','Linewidth',3);
grid on;
xlabel('ELA position x_0 (km from centre)');
ylabel('Ice sheet length (km)');
plot(x0,L,'ok','Linewidth',4,'Color','k','MarkerFaceColor','k', ...
    'MarkerSize',12);
grid on;
axis([Xmin Xmax Ymin Ymax]);
grid on;
plotval=strcat('x_0 = ',num2str(x0));
text(Xmin,Ymin-1800,plotval);
plotval=strcat('\beta = ',num2str(beta));
text(-520,Ymin-1800,plotval);
plotval=strcat('\epsilon = ',num2str(epsilon));
text(-250,Ymin-1800,plotval);
plotval=strcat('L = ',num2str(L),' km');
text(-500,Ymin-1200,plotval);
plotval=strcat('ELA = ',num2str(hf),' m a.s.l.');
text(-300,Ymin-1200,plotval);
text(-780,4700,'\bf \fontsize{12} Stable branch','BackgroundColor','w','Color','b');
text(-780,4400,'\bf \fontsize{12} Unstable branch','BackgroundColor','w','Color','r');
hold off;

end

