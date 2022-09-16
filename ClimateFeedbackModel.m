function ClimateFeedbackModel

%----------------------------------
% Frank Pattyn, 2016
% Laboratoire de Glaciologie, ULB
%----------------------------------


clear all;
close all

f = figure('Name','Climate Feedback Model','numbertitle','off', ...
    'Position', [200, 100, 900, 600]);

% Parameters
alb0=0.3;
albmax=0.6;
sigma=5.67e-8;
epsilon0=0.62;
So=1370;
Tice=-10;
Tnoice=0;
Tinit=15;
albice=alb0;


% Controls
Emin=50;
Emax=350;
Tmin=-50;
Tmax=60;
npts=100;
hysteresis=0;

T=linspace(Tmin,Tmax,npts);
Sfx=linspace(0.5,1.5,npts);
Epsx=linspace(0.4,1,npts);

%// Make initial plot
epsilon=epsilon0;
Sf=1;

Eout=EnergyOut(epsilon,sigma,T);
[Te,Ee]=eqvalue(Tinit,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma);
[Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So);
p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);

%// re-position the axes to make room for the sliders/buttons
set(gca, 'position', [0.1 0.25 0.8 0.7]);

%// initialize the sliders/buttons
h = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.05 0.05 0.25 0.05],...
    'min'     , 0.5,...               %// Make the A between 1...
    'max'     , 1.5,...              %// and 10, with initial value
    'value'   , Sf,...               %// as set above.
    'callback', @sliderSf);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar

                                    
h2=uicontrol(...
    'parent'  , f,...        
    'Style', 'popup',...
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'String', 'no ice|ice',...
    'Position', [0.6 0.05 0.08 0.05],...
    'Callback', @iceControl);
       
h3 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.32 0.05 0.25 0.05],...
    'min'     , 0.4,...               %// Make the A between 1...
    'max'     , 1,...              %// and 10, with initial value
    'value'   , epsilon,...               %// as set above.
    'callback', @sliderEpsilon);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar

h4 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'pushbutton',...    
    'String'  , 'INIT',...
    'position', [0.85 0.05 0.08 0.05],...
    'callback', @pushReset);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    

h5 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...    
    'String'  , 'Energy|Equilibrium Sf|Equilibrium E',...
    'position', [0.7 0.05 0.10 0.05],...
    'callback', @popupHyst);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar  
                                    
                                    
%// THE MAGIC INGREDIENT
%// ===========================

% hLstn = handle.listener(h,'ActionEvent',@sliderSf); %#ok<NASGU>
% hLstn = handle.listener(h2,'ActionEvent',@iceControl); %#ok<NASGU>
% hLstn = handle.listener(h3,'ActionEvent',@sliderEpsilon); %#ok<NASGU>
% hLstn = handle.listener(h4,'ActionEvent',@pushReset); %#ok<NASGU>
% hLstn = handle.listener(h5,'ActionEvent',@popupHyst); %#ok<NASGU>

%// (variable appears unused, but not assigning it to anything means that 
%// the listener is stored in the 'ans' variable. If "ans" is overwritten, 
%// the listener goes out of scope and is thus destroyed, and thus, it no 
%// longer works.

%// ===========================


%// The slider's callback:
%//    1) clears the old plot
%//    2) computes new values using the (continuously) updated slider values
%//    3) re-draw the plot and re-set the axes settings

function sliderSf(~,~)
    delete(p);
    Sf=get(h,'value');
    Eout=EnergyOut(epsilon,sigma,T);
    [Te,Ee]=eqvalue(Te,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma);
    [Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So);
    if hysteresis==0
        p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);
    else
        if hysteresis==1
            [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts);
            p=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts);
        else
            [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts);
            p=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts);
        end
    end
end

function sliderEpsilon(~,~)
    delete(p);
    epsilon=get(h3,'value');
    Eout=EnergyOut(epsilon,sigma,T);
    [Te,Ee]=eqvalue(Te,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma);
    [Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So);
    if hysteresis==0
        p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);
    else
        if hysteresis==1
            [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts);
            p=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts);
        else
            [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts);
            p=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts);
        end
    end
end

function iceControl(~,~)
    delete(p);
    iceval=get(h2,'value');
    Sf=get(h,'value');
    if iceval==1
        albice=alb0;
    else
        albice=albmax;
    end
    Eout=EnergyOut(epsilon,sigma,T);
    [Te,Ee]=eqvalue(Te,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma);
    [Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So);
    if hysteresis==0
        p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);
    else
        if hysteresis==1
            [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts);
            p=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts);
        else
            [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts);
            p=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts);
        end
    end
end

function pushReset(~,~)
    delete(p);
    epsilon=epsilon0;
    set(h3,'Value',epsilon0);
    Sf=1;
    Te=Tinit;
    set(h,'Value',1);
    Eout=EnergyOut(epsilon,sigma,T);
    [Te,Ee]=eqvalue(Te,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma);
    [Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So);
    if hysteresis==0
        p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);
    else
        if hysteresis==1
            [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts);
            p=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts);
        else
            [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts);
            p=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts);
        end
    end
end

function popupHyst(~,~)
    delete(p);
    graphtype=get(h5,'value');
    if graphtype==1
        hysteresis=0;
    else
        if graphtype==2
            hysteresis=1;
        else
            hysteresis=2;
        end
    end
    if hysteresis==0
        p=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon);
    else
        if hysteresis==1
            [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts);
            p=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts);
        else
            [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts);
            p=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts);
        end
    end
end


end

function [Ein,alb]=EnergyIn(npts,T,alb0,albice,Tice,Tnoice,Sf,So)

    alb=zeros(npts,1);
    for i=1:npts
        alb(i)=albedo(T(i),alb0,albice,Tice,Tnoice);
    end
    Ein = Sf*So*(1-alb)/4;

end

function [Eout]=EnergyOut(epsilon,sigma,T)

    Eout = epsilon*sigma*(T+273.15).^4;

end

function [p]=plotfigure(Sf,Te,Ee,Ein,Eout,T,Emin,Emax,Tmin,Tmax,epsilon)

    p = plot(T, Ein,'-r','Linewidth',3); hold on;
    plot(T,Eout,'-b','Linewidth',3);
    Tscale=(-90:1:90)';
    Tcolor=colormap(jet(length(Tscale)));
    plot(Te,Ee,'ok','MarkerFaceColor',Tcolor(Tscale==round(Te),:), ...
        'MarkerSize',25,'Linewidth',2);
    caxis([-90,90]);
%     colorbar;
    axis tight
    axis([Tmin Tmax Emin Emax])
    grid on;
    plotval=strcat('\epsilon = ',num2str(epsilon));
    text(-15,Emin-50,plotval);
    plotval=strcat('S_f = ',num2str(Sf));
    text(-50,Emin-50,plotval);
    plotval=strcat('T_{Earth} = ',num2str(Te),' °C');
    text(20,Emin-45,plotval);
    xlabel('Atmospheric temperature (°C)');
    ylabel('Energy (W m^{-2})');
    text(-48,330,'\bf \fontsize{12} Incoming solar radiation','BackgroundColor','w','Color','r');
    text(-48,315,'\bf \fontsize{12} Outgoing terrestrial radiation','BackgroundColor','w','Color','b');

    hold off;
end


function [T,E]=eqvalue(Tinit,Sf,So,alb0,albice,Tice,Tnoice,epsilon,sigma)

% determine equilibrium climate

T=(Sf*So*(1-albedo(Tinit,alb0,albice,Tice,Tnoice))/ ...
    (4*epsilon*sigma)).^0.25-273.15;
for i=1:10
    T=(Sf*So*(1-albedo(T,alb0,albice,Tice,Tnoice))/ ...
    (4*epsilon*sigma)).^0.25-273.15;
end
E=epsilon*sigma*(T+273.15).^4;

end

function alb=albedo(T,alb0,albice,Tice,Tnoice)

alb(T>Tnoice)=alb0;
alb(T<=Tnoice & T<Tice)=albice;
alb(T<=Tnoice & T>=Tice)=alb0+(albice-alb0)/(Tice-Tnoice)* ...
    (T(T<=Tnoice & T>=Tice)-Tnoice);
end

function [T1,T2]=HysteresisSf(Sfx,So,alb0,albice,epsilon,sigma,Tice,Tnoice,npts)

% Determine hysteresis curves
T1=zeros(npts,1)+80;
T2=zeros(npts,1)-80;
for j=1:5
    T1=(Sfx*So.*(1-albedo(T1,alb0,albice,Tice,Tnoice))/ ...
        (4*epsilon*sigma)).^0.25-273.15;
    T2=(Sfx*So.*(1-albedo(T2,alb0,albice,Tice,Tnoice))/ ...
        (4*epsilon*sigma)).^0.25-273.15;
end
if albice>alb0
    T1(T1<=Tice)=NaN;
    T2(T2>=Tice)=NaN;
end

end


function [p]=HysteresisSfplot(Sfx,T1,T2,Te,Sf,epsilon,npts)

% Makes hysteresis plot

p=plot(Sfx,T2,'b-','LineWidth',3); hold on;
plot(Sfx,T1,'r-','LineWidth',3);
B1=NaN(npts,1);
B2=NaN(npts,1);
for i=2:npts
    if isnan(T2(i)) && ~isnan(T2(i-1))
        B2(i-1)=T2(i-1);
        B2(i)=T1(i);
    end
    if ~isnan(T1(i)) && isnan(T1(i-1))
        B1(i)=T1(i);
        B1(i-1)=T2(i-1);
    end
end
plot(Sfx,B1,'r--','LineWidth',3);
plot(Sfx,B2,'b--','LineWidth',3);
plot(Sf,Te,'ok','MarkerFaceColor','k', ...
    'MarkerSize',15,'Linewidth',2);
grid on;
axis([0.5 1.5 -70 70]);
plotval=strcat('\epsilon = ',num2str(epsilon));
text(0.9,-95,plotval);
plotval=strcat('S_f = ',num2str(Sf));
text(0.6,-95,plotval);
plotval=strcat('T_{Earth} = ',num2str(Te),' °C');
text(1.2,-85,plotval);
xlabel('S_f');
ylabel('Atmospheric temperature (°C)');
text(0.52,60,'\bf \fontsize{12} Branch \alpha_{no ice}','BackgroundColor','w','Color','r');
text(0.52,50,'\bf \fontsize{12} Branch \alpha_{ice}','BackgroundColor','w','Color','b');

hold off;

end

function [T1,T2]=HysteresisEps(Epsx,So,alb0,albice,Sf,sigma,Tice,Tnoice,npts)

% Determine hysteresis curves
T1=zeros(npts,1)+80;
T2=zeros(npts,1)-80;
for j=1:5
    T1=(Sf*So*(1-albedo(T1,alb0,albice,Tice,Tnoice))./ ...
        (4*Epsx*sigma)).^0.25-273.15;
    T2=(Sf*So*(1-albedo(T2,alb0,albice,Tice,Tnoice))./ ...
        (4*Epsx*sigma)).^0.25-273.15;
end
if albice>alb0
    T1(T1<=Tice)=NaN;
    T2(T2>=Tice)=NaN;
end

end


function [p]=HysteresisEpsplot(Epsx,T1,T2,Te,Sf,epsilon,npts)

% Makes hysteresis plot

p=plot(Epsx,T2,'b-','LineWidth',3); hold on;
plot(Epsx,T1,'r-','LineWidth',3);
B1=NaN(npts,1);
B2=NaN(npts,1);
for i=2:npts
    if isnan(T1(i)) && ~isnan(T1(i-1))
        B1(i-1)=T1(i-1);
        B1(i)=T2(i);
    end
    if ~isnan(T2(i)) && isnan(T2(i-1))
        B2(i)=T2(i);
        B2(i-1)=T1(i-1);
    end
end
plot(Epsx,B1,'r--','LineWidth',3);
plot(Epsx,B2,'b--','LineWidth',3);
plot(epsilon,Te,'ok','MarkerFaceColor','k', ...
    'MarkerSize',15,'Linewidth',2);
grid on;
axis([0.4 1 -70 70]);
plotval=strcat('\epsilon = ',num2str(epsilon));
text(0.65,-95,plotval);
plotval=strcat('S_f = ',num2str(Sf));
text(0.45,-95,plotval);
plotval=strcat('T_{Earth} = ',num2str(Te),' °C');
text(0.85,-85,plotval);
xlabel('\epsilon');
ylabel('Atmospheric temperature (°C)');
text(0.87,60,'\bf \fontsize{12} Branch \alpha_{no ice}','BackgroundColor','w','Color','r');
text(0.87,50,'\bf \fontsize{12} Branch \alpha_{ice}','BackgroundColor','w','Color','b');

hold off;

end

