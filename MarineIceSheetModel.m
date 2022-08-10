function MarineIceSheetModel

% Schoof experiments steady-state

close all;
clear all;

f = figure('Name','Marine Ice Sheet Model','numbertitle','off', ...
    'Position', [200, 100, 900, 600]);

% Constants
secperyear=31556926;
a=0.3/secperyear;
A=1e-25;
n=3;
r=0.9; %ratio of ice to water density
rho_g=8820;
m=1/n;
C=7.624e6;
sl=0;
btype=0;
hysteresis=0;


% bed
xs=7.5e5;
b1=-720;
b2=900/xs;
a1=-729;
a2=2184.8/xs^2;
a3=-1031.72/xs^4;
a4=151.72/xs^6;

npts=500;
user_grid=linspace(0,1600e3,npts);
sealevel=linspace(150,-150,npts);
massb=linspace(0.05,0.6,npts)/secperyear;
grl0=NaN(1,length(sealevel));
grl1=grl0;
grl2=grl0;
grl3=grl0;
pos=10;  %initial guess of grounding line position (metres)


%// Make initial plot
[x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);

%// re-position the axes to make room for the sliders/buttons
set(gca, 'position', [0.1 0.25 0.8 0.7]);


%// initialize the sliders/buttons
h1 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.05 0.05 0.2 0.05],...
    'min'     , -150,...               %// Make the A between 1...
    'max'     , 150,...              %// and 10, with initial value
    'value'   , sl,...               %// as set above.
    'callback', @sliderSL);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h2 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'slider',...        
    'position', [0.28 0.05 0.2 0.05],...
    'min'     , 0.05,...               %// Make the A between 1...
    'max'     , 0.6,...              %// and 10, with initial value
    'value'   , a*secperyear,...               %// as set above.
    'callback', @sliderA);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
                                    
h3 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...        
    'position', [0.51 0.05 0.12 0.05],...
    'String'  , 'Linear|Overdeepened',...
    'callback', @popupLin);   %// This is called when using the arrows
                                    %// and/or when clicking the slider bar
                                    
h4 = uicontrol(...
    'parent'  , f,...        
    'units'   , 'normalized',...    %// so yo don't have to f*ck with pixels
    'style'   , 'popup',...    
    'String'  , 'Ice sheet|Equilibria|Phase space SL|Phase space a',...
    'position', [0.65 0.05 0.12 0.05],...
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

% hLstn = handle.listener(h1,'ActionEvent',@sliderSL); %#ok<NASGU>
% hLstn = handle.listener(h2,'ActionEvent',@sliderA); %#ok<NASGU>
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

function sliderSL(~,~)
    sl=get(h1,'value');
    [x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
    if hysteresis==0
        p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);
    else
        if hysteresis==1
            p=PlotEquilibria(x,user_grid,sl);
        else
            if hysteresis==2
                p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1);
            else
                p=PlotHysteresisA(x,sl,a,massb,grl2,grl3);
            end
        end
    end
end

function sliderA(~,~)
    a=get(h2,'value');
    a=a/secperyear;
    [x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
    if hysteresis==0
        p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);
    else
        if hysteresis==1
            p=PlotEquilibria(x,user_grid,sl);
        else
            if hysteresis==2
                p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1);
            else
                p=PlotHysteresisA(x,sl,a,massb,grl2,grl3);
            end
        end
    end
end

function popupLin(~,~)
    getLin=get(h3,'value');
    pos=10;
    if getLin==1
        btype=0;
    else
        btype=1;
    end
    [x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
    [sealevel,grl0,grl1]=HysteresisSL(user_grid,sealevel,grl0,grl1);
    [massb,grl2,grl3]=HysteresisA(user_grid,massb,grl2,grl3);
    if hysteresis==0
        p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);
    else
        if hysteresis==1
            p=PlotEquilibria(x,user_grid,sl);
        else
            if hysteresis==2
                p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1);
            else
                p=PlotHysteresisA(x,sl,a,massb,grl2,grl3);
            end
        end
    end
end

function pushInit(~,~)
    delete(p);
    set(h1,'Value',0);
    set(h2,'Value',0.3);
    sl=0;
    a=0.3/secperyear;
    pos=10;
    [x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
    [sealevel,grl0,grl1]=HysteresisSL(user_grid,sealevel,grl0,grl1);
    [massb,grl2,grl3]=HysteresisA(user_grid,massb,grl2,grl3);
    if hysteresis==0
        p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);
    else
        if hysteresis==1
            p=PlotEquilibria(x,user_grid,sl);
        else
            if hysteresis==2
                p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1);
            else
                p=PlotHysteresisA(x,sl,a,massb,grl2,grl3);
            end
        end
    end
end

function popupHyst(~,~)
    delete(p);
    graphtype=get(h4,'value');
    if graphtype==1
        hysteresis=0;
        [x,pos,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(pos,user_grid,sl);
        p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating);
    else
        if graphtype==2
            hysteresis=1;
            p=PlotEquilibria(x,user_grid,sl);
        else
            if graphtype==3
                hysteresis=2;
                [sealevel,grl0,grl1]=HysteresisSL(user_grid,sealevel,grl0,grl1);
                p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1);
            else
                hysteresis=3;
                [massb,grl2,grl3]=HysteresisA(user_grid,massb,grl2,grl3);
                p=PlotHysteresisA(x,sl,a,massb,grl2,grl3);
            end
        end
    end
end


function p=PlotEquilibria(x,user_grid,sl)
    
    b=bedheight(user_grid,sl);
    b(b<0)=NaN;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.));
    p=plot(user_grid/1e3,flux*secperyear/1e6,'r','Linewidth',3); hold on;
    plot(user_grid/1e3,a*user_grid*secperyear/1e6,'b','Linewidth',3);
    plot(x/1e3,a*x*secperyear/1e6,'o','Color','k','MarkerFaceColor','k', ...
            'MarkerSize',12);
    ylim([0 1]);
    xlim([0 1600]);
    grid on;
    xlabel('Grounding-line position (km)');
    ylabel('Ice flux (km^2 a^{-1})');

    plotval=strcat('SL = ',num2str(sl),' m');
    text(50,-0.18,plotval);
    plotval=strcat('a = ',num2str(a*secperyear),' m a^{-1}');
    text(500,-0.18,plotval);
    plotval=strcat('x_g = ',num2str(x*1e-3),' km');
    text(1200,-0.1,plotval);
    text(50,0.95,'\bf \fontsize{12} a.x','BackgroundColor','w','Color','b');
    text(50,0.88,'\bf \fontsize{12} Flux_{gl}','BackgroundColor','w','Color','r');
    hold off;

end

function p=PlotIceSheet(x,user_grid,sl,H_soln,X_soln,h_f,H_soln2,floating)
        
    p=plot(user_grid/1e3,-bedheight(user_grid,sl),'k','Linewidth',3); hold on;
    S_soln = H_soln - bedheight(X_soln,sl);
    plot(X_soln/1e3,S_soln,'b','Linewidth',3);
    plot([x/1e3; floating/1e3],(1-r)*[h_f; H_soln2],'b','Linewidth',3);
    plot([x/1e3; floating/1e3],-r*[h_f; H_soln2],'b','Linewidth',3);
    plot(x/1e3,-bedheight(x,sl),'o','Color','k','MarkerFaceColor','k', ...
            'MarkerSize',12);
    grid on;
    ylim([-3000 5500]);
    xlim([0 1600]);
    xlabel('Horizontal distance (km)');
    ylabel('Elevation (m above sea level)');

    plotval=strcat('SL = ',num2str(sl),' m');
    text(50,-4500,plotval);
    plotval=strcat('a = ',num2str(a*secperyear),' m a^{-1}');
    text(500,-4500,plotval);
    plotval=strcat('x_g = ',num2str(x*1e-3),' km');
    text(1200,-4000,plotval);
    text(1300,4900,'\bf \fontsize{12} Ice sheet','BackgroundColor','w','Color','b');
    text(1300,4500,'\bf \fontsize{12} Bedrock','BackgroundColor','w','Color','k');
    hold off;

end


function p=PlotHysteresisSL(x,sl,sealevel,grl0,grl1)
    
    p=plot(sealevel,grl1/1e3,'r','Linewidth',3); hold on;
    plot(sealevel,grl0/1e3,'b','Linewidth',3);
    plot(sl,x/1e3,'o','Color','k','MarkerFaceColor','k', ...
            'MarkerSize',12);
    grid on;
    ylim([600 1500]);
    xlabel('Sea level (m)');
    ylabel('Grounding-line position (km)');

    plotval=strcat('SL = ',num2str(sl),' m');
    text(-140,450,plotval);
    plotval=strcat('a = ',num2str(a*secperyear),' m a^{-1}');
    text(-60,450,plotval);
    plotval=strcat('x_g = ',num2str(x*1e-3),' km');
    text(80,500,plotval);
    text(-140,700,'\bf \fontsize{12} Advance branch','BackgroundColor','w','Color','b');
    text(-140,650,'\bf \fontsize{12} Return branch','BackgroundColor','w','Color','r');

    hold off;
end


function p=PlotHysteresisA(x,sl,a,massb,grl2,grl3)
    
    p=plot(massb*secperyear,grl3/1e3,'r','Linewidth',3); hold on;
    plot(massb*secperyear,grl2/1e3,'b','Linewidth',3);
    plot(a*secperyear,x/1e3,'o','Color','k','MarkerFaceColor','k', ...
            'MarkerSize',12);
    grid on;
    xlim([0.05 0.6]);
    ylim([600 1500]);
    xlabel('Accumulation rate (m a^{-1})');
    ylabel('Grounding-line position (km)');

    plotval=strcat('SL = ',num2str(sl),' m');
    text(0.07,450,plotval);
    plotval=strcat('a = ',num2str(a*secperyear),' m a^{-1}');
    text(0.22,450,plotval);
    plotval=strcat('x_g = ',num2str(x*1e-3),' km');
    text(0.47,500,plotval);
    text(0.47,700,'\bf \fontsize{12} Advance branch','BackgroundColor','w','Color','b');
    text(0.47,650,'\bf \fontsize{12} Return branch','BackgroundColor','w','Color','r');

    hold off;
end


function [sealevel,grl0,grl1]=HysteresisSL(user_grid,sealevel,grl0,grl1)
    
i=10;
for j=1:length(sealevel)
    flag=0;
    b=bedheight(user_grid,sealevel(j));
    b(b<0)=0;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.))-a*user_grid;
    while flag==0
        if flux(i)<0
            i=i+1;
        else
            i=i-1;
        end
        if flux(i)<=0 && flux(i+1)>0
            flag=1;
        end
    end
    grl0(j)=user_grid(i);
end
for j=length(sealevel):-1:1
    flag=0;
    b=bedheight(user_grid,sealevel(j));
    b(b<0)=0;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.))-a*user_grid;
    while flag==0
        if flux(i)<0
            i=i+1;
        else
            i=i-1;
        end
        if flux(i)<=0 && flux(i+1)>0
            flag=1;
        end
    end
    grl1(j)=user_grid(i);
end
    
end

function [massb,grl2,grl3]=HysteresisA(user_grid,massb,grl2,grl3)
    
i=10;
for j=1:length(massb)
    flag=0;
    b=bedheight(user_grid,sl);
    b(b<0)=0;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.))-massb(j)*user_grid;
    while flag==0
        if flux(i)<0
            i=i+1;
        else
            i=i-1;
        end
        if flux(i)<=0 && flux(i+1)>0
            flag=1;
        end
    end
    grl2(j)=user_grid(i);
end
for j=length(massb):-1:1
    flag=0;
    b=bedheight(user_grid,sl);
    b(b<0)=0;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.))-massb(j)*user_grid;
    while flag==0
        if flux(i)<0
            i=i+1;
        else
            i=i-1;
        end
        if flux(i)<=0 && flux(i+1)>0
            flag=1;
        end
    end
    grl3(j)=user_grid(i);
end
    
end


function [x,i,X_soln,h_f,H_soln,H_soln2,floating]=IceSheetModel(i,user_grid,sl)

    flag=0;
    b=bedheight(user_grid,sl);
    b(b<0)=0;
    flux=((A*(rho_g)^(n+1.)*(1-r)^n)/(4.^n*C))^(1./(m+1.))* ...
        (b/r).^((m+n+3.)/(m+1.))-a*user_grid;
    while flag==0
        if flux(i)<0
            i=i+1;
        else
            i=i-1;
        end
        if flux(i)<=0 && flux(i+1)>0
            flag=1;
        end
    end
    x=user_grid(i);

    grounded = user_grid(user_grid < x).';
    x_grid = [x; grounded(length(grounded):-1:1)];
    h_f = bedheight(x,sl)/r;
    options = odeset('AbsTol',1e-6,'RelTol',1e-6); %odeset('AbsTol',f/1e-3);
    [X_soln,H_soln] = ode45(@SMsurface,x_grid,[h_f],options);
    floating = user_grid(user_grid > x).';
    q_0 = a*x;
    H_soln2 = h_f*(q_0 + a*(floating-x))./(q_0^(n+1) + h_f^(n+1) ...
        *((1-r)*rho_g/4)^n*A*((q_0 + a*(floating-x)).^(n+1) ...
        -q_0^(n+1))/a).^(1/(n+1));
    H_soln2(floating>x+150e3)=NaN;

end


function z = bedheight(x,sl)
%DEPTH OF BED BELOW SEA LEVEL, i.e. gives b(x) in Schoof 2007.
%NOTE the sign convention, SMcold_bedheight is positive if the bed is below
%sea level, negative if above sea level.

if btype==1
    z=sl+a1+a2*x.^2 +a3*x.^4+a4*x.^6;
else
    z=sl+b1+b2*x;
end

end


function z = bedslope(x)
%FIRST DERIVATIVE OF DEPTH OF BED BELOW SEA LEVEL; must agree with
%SMcold_bedheight.

if btype==1
    z=2*a2*x+4*a3*x.^3+6*a4*x.^5;
else
    z=b2+x*0;
end

end


function z = SMsurface(x,h)

b_x = bedslope(x);
s = a*x;

z = b_x - (C/rho_g)*s.^m./h.^(m+1);

end

end




