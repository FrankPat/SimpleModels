%-----------------------------------------------
% Simple Climate and Sea Level Model (SCM-SLRM)
%-----------------------------------------------
%
% Simple radiative model based that predicts global mean temperature as a
% function of CO2 equivalent atmosphere concentrations. The latter is the
% basic input for the model. CO2-eq. is related to radiative forcing as:
%
%   dQ = 5.35 ln([CO2]/[CO2_o])
%
% where [CO2_o] is the preindustrial level (280 ppmv). The black body
% temperature is then:
%
%   dTs0 = -dQ/L0, where L0 = -4 \sigma Te^3
%
% where \sigma is the Stefan-Boltzmann constant and Te the temperature at
% the top of the troposphere (255 K). Atmospheric temperartures are then
% related to the black body temperature via a feedback factor Ff:
%
%   dTs = Ff . dTs0
%
% where Ff = La0/Laf and Laf is the sum of all feedbacks. Ff ~ 1.74. Ff can
% also be determined analytically by correlating dTs0 directly to observed
% temperatures (OptimizeFeedback = 1).
%
% Future scenarios rely on the relationship between carbon emissions and
% CO2-eq. concentrations in the atmosphere, taken from the observed changes
% in the last 50 years. The relation is automatically calculated.
%
% Datasets used are:
% (1) observed CO2-eq. atmosphere concentration, taken into
% account CO2, CH4, N20 emissions as well as negative forcing effects from
% aerosols.
% (2) Historical CO2 emissions, including land-use changes (Gt CO2)
% (3) Historical temperature changes (for optimization of Ff)
% (4) Future carbon emissions scenarios (SSPs)
%
%-----------------------------------------------
% DATA
%-----------------------------------------------
%
% Mean global surface temperatures (yearly)
% This file contains a brief summary of the changes in Earth's global average
% surface temperature estimated by combining the Berkeley Earth land-surface
% temperature field with a reinterpolated version of the HadSST4 ocean temperature 
% field.  
% The current citation for this dataset is: 
%    Rohde, R. A. and Hausfather, Z.: The Berkeley Earth Land/Ocean Temperature
%    Record, Earth Syst. Sci. Data, 12, 3469-3479, 
%    https://doi.org/10.5194/essd-12-3469-2020, 2020.
% The dataset differs slightly from the dataset as described in the citation as 
% HadSST3 has been replaced with the newer HadSST4, and associated interpolation 
% parameters have been refit accordingly.  No other changes in methods were needed 
% when moving to the new version of HadSST. 
% http://berkeleyearth.lbl.gov/auto/Global/Land_and_Ocean_summary.txt
%
% Observed trends in total greenhouse gas concentration levels between 1860 and 2019,
% considering all greenhouse gases and other forcing agents (including aerosols)
% https://www.eea.europa.eu/ims/atmospheric-greenhouse-gas-concentrations
%
% Global CO2 emissions from fossil fuels
% https://ourworldindata.org/co2-emissions
%
% SSP emission scenarios:
% https://www.ipcc.ch/data/
% https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=10
%
% To convert from gigatonnes carbon to gigatonnes of carbon dioxide, you
% simply multiply 44 over 12. In other words, 1 gigatonne of carbon
% equals 3.67 gigatonnes of carbon dioxide


function SCM_SLRM(scenario)

close all;

% Control parameters
startyeartmp=1850;
currentyear=2019;
endyear=2100; % En year of simulation
OptimizeFeedback=0; % Determine the feedback factor Ff from the relation
                    % between dTs0 and observed temperatures
LandUseChange=1; % Incorporate land use change in carbon emissions
NoFeedback=0;
OceanHeat=1;

% Sea level control parameters
sealevel=1;
startyearsea=1880;
ed = 15; % embedding dimension (years) of SSAtrend smoother (15 was used in the paper)
iterate = 'y'; % iterate for lambda. No iteration -> estimate a only
bin=15;
chao = 'y'; % use Chao reservoir corr

% Scenario
scenario_txt={'Observed','SSP1-19','SSP1-26','SSP2-45', ...
    'SSP3-70','SSP5-85','Constant','T2xCO2'};
if nargin<1
    fprintf('Future scenario:\n');
    for i=1:8
        fprintf('  (%d) %s\n',i-1,char(scenario_txt(i)));
    end
    scenario=input('Choose scenario: ');
end

% constants
boltzmann=5.67e-8; % Stefan-Boltzmann constant
albedo=0.3; % Planetary albedo
solar=1370; % Solar constant
CO2o=280; % Pre-industrial CO2 concentration
Cs=8.36e8/31536000; % characteristic time scale thermal inertia

% feedbacks
Te=(solar*(1-albedo)/(4*boltzmann))^0.25;
La0=-4*boltzmann*Te^3; % Black body radiation feedback
LaW=1.6; % 1.6: water vapour feedback
LaL=-0.6; % -0.6: Lapse rate feedback
LaC=0.3; %-0.4 - +1: Cloud feedback (uncertain)
LaA=0.3; % 0.3 % albedo feedback
LaUnc=0.5; % uncertainty on feedback parameters

if OceanHeat==1 % data from IPCC AR6
    LaW=1.6; % 1.6: water vapour feedback
    LaL=-0.3; % -0.6: longwave radiation feedback
    LaC=0.42; %-0.4 - +1: Cloud feedback (uncertain)
    LaA=0.35; % 0.3 % albedo feedback
    LaUnc=0.65; % uncertainty on feedback parameters
end

Lf=La0+LaW+LaL+LaC+LaA;
if OceanHeat==1
    tau=-Cs/Lf;
else
    tau=0;
end

if NoFeedback==1
    OptimizeFeedback=0;
    LaW=0;
    LaL=0;
    LaC=0;
    LaA=0;
end

% Read input files

data=load('GreenhouseGases.txt');
year=(startyeartmp:1:currentyear)';
CO2=interp1(data(:,1),data(:,2),year);

data=load('BerkeleyTemps2021.txt');
tmp=interp1(data(:,1),data(:,6),year);
% data=load('HadCrutTemps.txt');
% tmp=interp1(data(:,1),data(:,2),year);
tmp=tmp-mean(tmp(1:70));

data=load('Emissions.txt');
if LandUseChange==1
    Em=interp1(data(:,1),data(:,2)/1e9,year);
else
    Em=interp1(data(:,1),data(:,4)/1e9,year);
end
Emtot=Em;
for i=2:length(year)
    Emtot(i)=Emtot(i-1)+Em(i);
end

% Fit between carbon emissions and CO2 concentrations
P=polyfit(Emtot(120:end),CO2(120:end),1);

data=load('SSPscenarios.txt');
year2=(currentyear:1:endyear)';
data(1,1)=currentyear;
if endyear>2100
    data(:,end+1)=data(:,end);
    data(1,end)=endyear;
end
data(2:end,1)=Em(end)*1e3;
SSP=zeros(length(year2),5);
for i=1:5
    SSP(:,i)=interp1(data(1,:),data(i+1,:)/1e3,year2);
end

if scenario>=7
    CO2(1:120)=CO2o;
    for i=121:170
        CO2(i)=CO2(i-1)*1.01;
    end
end

% Radiative forcing
dQ=CO2dQ(CO2,CO2o);
dTs0=-dQ/La0;
ft1=fittype({'x'});
if OptimizeFeedback==1
    Q=fit(dTs0,tmp,ft1);
    yfit=feval(Q,dTs0);
    Ff=feval(Q,1);
    FfUp=Ff+La0/(Ff*La0+LaUnc);
    FfLo=Ff-La0/(Ff*La0+LaUnc);
else
    Ff=1/(1+LaL/La0+LaW/La0+LaC/La0+LaA/La0);
    FfUp=1/(1+LaL/La0+LaW/La0+LaC/La0+LaA/La0+LaUnc/La0);
    FfLo=1/(1+LaL/La0+LaW/La0+LaC/La0+LaA/La0-LaUnc/La0);
    yfit=dTs0*Ff;
end
dTs=Ff*dTs0;
dTsUp=FfUp*dTs0;
dTsLo=FfLo*dTs0;
if tau>0
    for i=2:length(year)
        dTs(i)=(dTs(i-1)+dQ(i)/tau)/(1-La0/(tau*Ff));
        dTsUp(i)=(dTsUp(i-1)+dQ(i)/tau)/(1-La0/(tau*FfUp));
        dTsLo(i)=(dTsLo(i-1)+dQ(i)/tau)/(1-La0/(tau*FfLo));
    end
end

% Future simulations

switch scenario
    case 0
        year2=(currentyear:1:currentyear+1)';
        FutEm=linspace(Em(end),Em(end),length(year2))';
    case {1,2,3,4,5}
        FutEm=SSP(:,scenario);
    case {6,7}
        FutEm=linspace(Em(end),Em(end),length(year2))';
end

FutEmtot=zeros(size(FutEm))+Emtot(end);
for i=2:length(year2)
    FutEmtot(i)=FutEmtot(i-1)+FutEm(i);
end

FutCO2=FutEmtot*P(1)+P(2);
if scenario==7
    FutCO2(1)=CO2(end);
    for i=2:length(FutCO2)
       FutCO2(i)=min(FutCO2(i-1)*1.01,2*CO2o);
    end
end
FutdQ=CO2dQ(FutCO2,CO2o);
FutdTs=-Ff*FutdQ/La0;
FutdTsUp=-FfUp*FutdQ/La0;
FutdTsLo=-FfLo*FutdQ/La0;
if tau>0
    FutdTs(1)=dTs(end);
    FutdTsUp(1)=dTsUp(end);
    FutdTsLo(1)=dTsLo(end);
    for i=2:length(year2)
        FutdTs(i)=(FutdTs(i-1)+FutdQ(i)/tau)/(1-La0/(tau*Ff));
        FutdTsUp(i)=(FutdTsUp(i-1)+FutdQ(i)/tau)/(1-La0/(tau*FfUp));
        FutdTsLo(i)=(FutdTsLo(i-1)+FutdQ(i)/tau)/(1-La0/(tau*FfLo));
    end
end

% Figures

scrsz = get(groot,'ScreenSize');
figure('Position',[30 scrsz(4)/3 scrsz(3)/1.6 scrsz(4)/1.8]);
subplot(2,2,1);
p1=plot(Emtot(120:end),CO2(120:end),'or'); hold on;
set(p1,'markersize',5,'markerfacecolor','r');
plot(Emtot(120:end),Emtot(120:end)*P(1)+P(2),'--','linewidth',2);
txt=['CO_2 = ' num2str(P(1)) ' C_c + ' num2str(P(2)) ' ppmv'];
text(Emtot(120)+100,450,txt,'FontSize',10);
grid on;
xlabel('Cumulative C emissions (Gt CO_2)');
ylabel('Atmospheric CO_2-eq. (ppmv)');

subplot(2,2,2);
hold on;
xh1=(0:0.1:max(dTs0))';
xh=[xh1;flipud(xh1)];
yh=[xh1*FfUp;flipud(xh1*FfLo)];
fill(xh,yh,[.8 .8 1],'marker','.','markersize',0.01,'edgecolor',[.8 .8 1]);
plot(dTs0,tmp,'o');
plot(dTs0,yfit,'k','linewidth',2);
grid on;
txt=['F_f = ' num2str(Ff) ' '];
text(0.1,1,txt,'FontSize',10);
xlabel('\Delta T_{s,0}');
ylabel('Observed \Delta T_s');

subplot(2,2,3);
R=fit(Emtot,tmp,ft1);
yfit=feval(R,Emtot);
plot(Emtot,tmp,'o'); hold on;
plot(Emtot,dTs,'r+');
plot(FutEmtot,FutdTs,'r+');
plot(Emtot,yfit,'k','linewidth',2);
grid on;
txt=['1000 Gt = ' num2str(feval(R,1)*1e3) ' °C'];
text(600,-0.2,txt,'FontSize',10);
xlabel('Cumulative C emissions (Gt CO_2)');
ylabel('\Delta T_s');
legend('Observed','Simulated','location','northwest');

subplot(2,2,4);
plot(tmp,dTs,'o'); hold on;
grid on;
xx=[-0.5 1.5];
plot(xx,xx,'linewidth',2);
xlabel('Observed \Delta T_s');
ylabel('Simulated \Delta T_s');


scrsz = get(groot,'ScreenSize');
figure('Position',[30 scrsz(4)/3 scrsz(3)/1.6 scrsz(4)/1.8]);
subplot(2,2,1);
plot(year,Em,'linewidth',2); hold on;
plot(year2,FutEm,'linewidth',2);
plot(year(end),Em(end),'ko','MarkerFaceColor','k');
grid on;
xlabel('Year');
ylabel('Carbon emissions (Gt CO_2 a^{-1})');
legend('Observed','Future','location','northwest');
text(1870,20,char(scenario_txt(scenario+1)));

subplot(2,2,3);
plot(year,Emtot,'linewidth',2); hold on;
plot(year2,FutEmtot,'linewidth',2);
plot(year(end),Emtot(end),'ko','MarkerFaceColor','k');
grid on;
xlabel('Year');
ylabel('Cumulative C emissions (Gt CO_2)');
legend('Observed','Future','location','northwest');

subplot(2,2,2);
plot(year,CO2,'linewidth',2); hold on;
plot(year2,FutCO2,'linewidth',2);
plot(year(end),CO2(end),'ko','MarkerFaceColor','k');
grid on;
xlabel('Year');
ylabel('Atmospheric CO_2-eq. (ppmv)');
legend('Observed','Future','location','northwest');

subplot(2,2,4);
hold on; box on;
xh=[year; flipud(year)];
yh=[dTsUp; flipud(dTsLo)];
fill(xh,yh,[.8 .8 1],'marker','.','markersize',0.01,'edgecolor',[.8 .8 1]);
xh=[year2; flipud(year2)];
yh=[FutdTsUp; flipud(FutdTsLo)];
fill(xh,yh,[.8 .8 1],'marker','.','markersize',0.01,'edgecolor',[.8 .8 1]);
plot(year,tmp,'b','linewidth',2);
plot(year,dTs,'r','linewidth',2);
plot(year2,FutdTs,'r','linewidth',2);
plot(year(end),dTs(end),'ko','MarkerFaceColor','k');
grid on;
xlabel('Year');
ylabel('Temperature change (°C)');
txt=['\Delta T_s (' num2str(year2(end)) ') = ' num2str(FutdTs(end)) ' °C'];
text(1950,-0.25,txt,'FontSize',10);


%---------------------------------------------
% sea level rise
%---------------------------------------------
% This code produces results of Vermeer and Rahmstorf. 
% "Global Sea Level Linked to Global Temperature", PNAS 2009
% For a description of what it does, see that paper.
%
% Code written by Stefan Rahmstorf and Martin Vermeer
%
% Please note that this is a scientific code to be used by scientists
% who know what they are doing; it is not "fool-proof" for the general 
% user in the sense that all combinations of parameter choices, options etc. 
% have been thoroughly tested.
% Please report any problems or errors to the authors.

%------------------------------
% Parameter and options setting
%------------------------------
% The code uses the SSAtrend of Moore et al. 2005 for smoothing the data,
% which can be downloaded at http://www.glaciology.net/software/ssatrend-m
% (As a simple matlab alternative, fitting a polynomial is available as
% an option below. This is not recommended; it gives only a poor fit, 
% although it has little effect on the future projections.) 
% Reference: J. C. Moore, A. Grinsted, S. Jevrejeva, Eos 86, 226 (2005).

% Data: GMSL dataset at CSIRO:
% https://research.csiro.au/slrwavescoast/sea-level/measurements-and-data/sea-level-data/

% lambda: coefficient for linear combination of T and dT/dt (MV)
% dH/dt = a T + b dT/dt = a (T + lambda dT/dt), so lambda = b/a

if sealevel==1   
    lambda = 0;

    years=(startyearsea:1:currentyear)';

    data=load('CSIRO_Recons_gmsl_yr_2019.txt');
    SLR=interp1(data(:,1),data(:,2)-data(1,2),years);
    data=load('BerkeleyTemps.txt');
    tmp=interp1(data(:,1),data(:,6),years);

    % data=load('SeaLevelRise.txt');
    % SLR=interp1(data(:,1),data(:,2)-data(1,2),years);
    % data=load('HadCrutTemps.txt');
    % tmp=interp1(data(:,1),data(:,2),years);

    % data=load('PNAS/church_13221.txt');
    % SLR=interp1(data(:,1),data(:,2)-data(1,2),years);
    % data=load('PNAS/giss_landocean.txt');
    % tmp=interp1(data(:,1),data(:,14)*0.01,years);

    tmp=tmp-mean(tmp(1:70));

    % Apply Chao et al (2008) reservoir correction:
    if chao == 'y'
       SLR=SLR+1.65+(3.7/pi)*atan2(years-1978,13);
    end

    % apply smoothing to time series
%     smoothsealevel = ssatrend([years SLR],ed);
%     smoothtemp = ssatrend([years tmp],ed);
    smoothsealevel = smooth(SLR,ed,'moving');
    smoothtemp = smooth(tmp,ed,'moving');

    rateofrise=diff(smoothsealevel); % time derivative of smoothed sea level
    rateoftemp=diff(smoothtemp); % time derivative of smoothed temperature
    rateofrise(end+1)=rateofrise(end);
    rateoftemp(end+1)=rateoftemp(end);

    dataspan = 1:floor((length(SLR)-1)/bin); % for any binning
    mtemp=zeros(max(dataspan),1);
    mrtmp=zeros(max(dataspan),1);
    mrate=zeros(max(dataspan),1);
    temp=flip(smoothtemp);
    rtmp=flip(rateoftemp);
    rate=flip(rateofrise);
    for i=dataspan
        mtemp(i)=mean(temp((i-1)*bin+1:i*bin));
        mrtmp(i)=mean(rtmp((i-1)*bin+1:i*bin));
        mrate(i)=mean(rate((i-1)*bin+1:i*bin));
    end

    % loop for optimal lambda:
    rold = 0;
    [r,p] = corrcoef([mtemp(dataspan) mrate(dataspan)]);
    while r(1,2) > rold && iterate == 'y'
        rold = r(1,2);
        lambda = lambda - 0.01;
        % construct a linear combination of T and dT/dt (MV):
        mtemp2 = mtemp + lambda * mrtmp;
        % correlation of these binned data:
        [r,p] = corrcoef([mtemp2(dataspan) mrate(dataspan)]);
    end
    if iterate == 'y'
        mtemp = mtemp2;
    end

    disp([10, 'Pearson correlation coefficient r: ']); disp(r(1,2));  % correlation coefficient

    % linear fit to these data:
    [po, s] = polyfit(mtemp(dataspan),mrate(dataspan),1); % s is used with "polyval" to compute error bounds
    [fitrate,ferror] = polyval(po,mtemp(dataspan),s);
    fitrate = polyval(po,[min(tmp) max(tmp)],s);

    slope = po(1); % this is the a parameter
    tbase = -po(2)/po(1); % this is the derived pre-industrial base temperature
    b = slope * lambda; % this is the b parameter
    fprintf(1, 'T0 = %6.3f\n a = %6.3f\n b = %6.3f\n lambda = %6.3f\n', tbase, slope, b, lambda);

    % compute sea level hindcast from the fit found above 
    [hindrate, herror] = polyval(po,smoothtemp+lambda*rateoftemp,s); %hindcast rate of sea level rise
    hindlevel = cumsum(hindrate); %integrate up the rate

    % compute future sea level
    fac = 1/0.675; %convert herror to s.d.
    Futrateoftemp=diff(FutdTs); % time derivative of smoothed temperature
    Futrateoftemp(end+1)=Futrateoftemp(end);
    [Futrate, perror] = polyval(po,FutdTs+lambda*Futrateoftemp,s); %hindcast rate of sea level rise
    Futlevel = cumsum(Futrate)+hindlevel(end-1); %integrate up the rate
    predmax=cumsum(Futrate+fac*perror)+hindlevel(end-1);
    predmin=cumsum(Futrate-fac*perror)+hindlevel(end-1);

    % compute future sea level with temperature error
    FutrateUp=diff(FutdTsUp); % time derivative of smoothed temperature
    FutrateUp(end+1)=FutrateUp(end);
    [FutUp, perror] = polyval(po,FutdTsUp+lambda*FutrateUp,s); %hindcast rate of sea level rise
    FutlevelUp = cumsum(FutUp+fac*perror)+hindlevel(end-1); %integrate up the rate
    FutrateLo=diff(FutdTsLo); % time derivative of smoothed temperature
    FutrateLo(end+1)=FutrateLo(end);
    [FutLo, perror] = polyval(po,FutdTsLo+lambda*FutrateLo,s); %hindcast rate of sea level rise
    FutlevelLo = cumsum(FutLo-fac*perror)+hindlevel(end-1); %integrate up the rate


    % Figures

    scrsz = get(groot,'ScreenSize');
    figure('Position',[30 scrsz(4)/3 scrsz(3)/1.6 scrsz(4)/1.8]);
    % scatter plot of actual vs modelled rate of sea level rise

    subplot(2,2,1);
    hold on; box on;
    p1=plot(mtemp,mrate,'or');
    set(p1,'markersize',5,'markerfacecolor','r');
    plot([min(tmp) max(tmp)],fitrate,'--');
    grid on;
    ylabel('Rate of SL Change (mm a^{-1})')
    xlabel('\Delta T_s (K)')

    % plot of time series of actual and modelled rate of rise
    subplot(2,2,3);
    hold on; box on;
    %plot error band of modelled rate
    xh=[years; flipud(years)];
    % re-smoothing to take out diff artefacts:
    hindrate = ssatrend([years hindrate],ed);
    yh=[(hindrate+fac*herror); flipud((hindrate-fac*herror))];
    % This is the light blue error band
    if iterate == 'y'
        fill(xh,yh,[.8 .8 1],'marker','.','markersize',0.01,'edgecolor',[.8 .8 1]);
    else
        fill(xh,yh,[.8 .8 .8],'marker','.','markersize',0.01,'edgecolor',[.8 .8 .8]);
    end
    smoothrate = ssatrend([years rateofrise],ed);
    plot(years,hindrate,'b','linewidth',1.5);
    plot(years,smoothrate,'r','linewidth',2);
    ylabel('Rate of Change (mm a^{-1})');
    xlabel('Year');
    grid on;

    subplot(1,2,2);
    hold on; box on;
    SLR0=hindlevel(years==1990);
    %plot range as grey band
    xs=[year2; flipud(year2)];
    ys=[FutlevelUp; flipud(FutlevelLo)];
    fill(xs,ys-SLR0,[.8 .8 1],'marker','.','markersize',0.01,'edgecolor',[.8 .8 .8]);
    xs=[year2; flipud(year2)];
    ys=[predmax; flipud(predmin)];
    fill(xs,ys-SLR0,[.8 .8 .8],'marker','.','markersize',0.01,'edgecolor',[.8 .8 .8]);
    plot(years,SLR-SLR0,'b','linewidth',1.5);
    plot(years,hindlevel-SLR0,'r','linewidth',1.5);
    plot(year2,Futlevel-SLR0,'r','linewidth',1.5);
    plot(xlim,[0,0],'linewidth',.5,'color',[.5 .5 .5]);
    plot([1990 1990],ylim,'linewidth',.5,'color',[.5 .5 .5]);
    grid on;
    text(1900,50,char(scenario_txt(scenario+1)));
    ylabel('Sea level rise since 1990 (mm)')
    xlabel('Year')
end

save toto;

end



%---------------------------------------------

function dQ=CO2dQ(CO2,CO2o)
    dQ=5.35*log(CO2/CO2o);
end


function [R,V,sE_R,trendidx,tfilt]=ssatrend(x,M,sE_x,varargin)
% Finds the ssa non-linear trend from a time series.
% 
%
% USAGE:
%     [R,V,sE_R,trendidx,tfilt]=ssatrend(x,M[,sE_x,parameterlist])
%
% INPUT:
%   x: the time series
%   M: embedding dimension (=maximum lag), (or a symmetric filter in the time domain, it will not be using ssa then)
%   SE_x: standard error of a measurement in x (assumend to be independent)
%         - can also be replaced by a function handle to a function that
%           generates noise-surrogates of the same length as x. E.g.
%           E.g: @()randn(size(x))*.6 is equivalent to .6
%
% PARAMETERS:
%   SSAMethod: see help on ssa
%   ExtrapMethod: MinimumRoughness(=Default) or MinimumSlope
%
% OUTPUT:
%   R: the nonlinear trend (a specific reconstructed component where
%      x has succesively padded with the linear trend of the last M points. For
%      details check source).
%      note that the mean(x) has been added to the RC.
%      (so that it is easy to plot on the same plot as the original
%      timeseries)
%   V: the variance explained by R (as a fraction of the total variance)
%   SE_R: standard error of the R (how much of sE_x is remaining in R).
%   trendidx: the index of the rc corresponding to the trend
%   tfilt: symmetric reconstructive SSA filter. (what is convoluted on x)
%
% if no output is specified then the function will make a plot where 2
% standard errors (2*sE_R) is shaded grey.
%
% (Requires ssa.m by Eric Breitenberger and formatts.m from my WTC toolbox)
%
% Note: the ssatrend filter very often results in virtually the same output
% as a simple triangular FIR filter of length 2*M-1.
%
% (c) Aslak Grinsted 2009
% 2009: added boostrapping estimates of padding errors. - please note the
% algorithm at present assumes white noise uncertainties of the trend.
%

Args=struct('SSAMethod','unbiased',...
    'ExtrapMethod','minimumroughness');
Args=parseArgs(varargin,Args);



if nargin<3, sE_x=[]; end

if nargin==0
    warning('No data specified, using debug data...')
    x=cumsum(randn(300,1)*.2)+randn(300,1);
    M=100;
    sE_x=1;
    Args.ExtrapMethod='mr';
end

x=formatts(x);
t=x(:,1);
x=x(:,2);
meanx=mean(x);
if ~isa(sE_x,'function_handle') %then make it into a function-handle
    sE_x=sE_x(:);
    %     if length(sE_x)==1
    %         sE_x=ones(size(x))*sE_x;
    %     end
    
    if length(sE_x)==1
        sE_x=@()randn(size(x))*sE_x;
    else
        sE_x=@()r;
    end
end

noerrorbar=1;
% noerrorbar=isempty(sE_x);
if noerrorbar
    sE_x=zeros(size(x,1),1);
end


if mean(sE_x())>std(x)*1.01
    warning('sE_x>std(x)!!!')
end

n=length(x);

if numel(M)==1
    
    [E,V,A,R,C]=ssa(x,M,Args.SSAMethod);
    
    
    
    %find trend idx:
    [mx,trendidx]=max(abs(sum(E))./sum(abs(E))); %trend would be non oscillatory and therefore have a large sum in the index
    
    
    %throw away everything but the trend:
    E=E(:,trendidx);
    V=V(trendidx)/sum(V); %convert to percent
    
    if V<.1
        warning('ssa trend accounts for less than 10% of the variance.')
    end
    R=R(:,trendidx)+meanx;
    tfilt=conv(E,flipud(E))/M;
else
    tfilt=M(:);
    V=nan;trendidx=nan;
    M=ceil(length(tfilt)/2); %used for padding.
end
%tfilt=tfilt/sum(tfilt); %NORMALIZE!!!!!!! (NO!)



%succesively pad the series with values from the trend and reconstruct
%again to get better end estimates.
% paddedX=x;
% for ii=1:M
%     leftR=mean(paddedX(1:M));
%     rightR=mean(paddedX(end+1-(1:M)));
%     paddedX=[leftR;paddedX;rightR];
% end
mp=M-1;
switch Args.ExtrapMethod
    case {'minimumroughness','minrough','mr'}
        idx=(1:mp)';
        
        %pre calculate all mp-windowed linear fits for the bootstrap. The
        %prediction residuals will be used to calculate the uncertainty due
        %to padding.
        pt=[idx ones(mp,1)];
        ptleft=[idx-mp ones(mp,1)];
        ptright=[idx+mp ones(mp,1)];
        residleft=nan(mp,n-mp*2);
        residright=nan(mp,n-mp*2);
        for ii=1:n-mp
            %pp(:,ii)=pt\x(idx+ii-1);
            pp=pt\x(idx+ii-1);
            if ii>mp
                residleft(:,ii-mp)=ptleft*pp-x(idx+ii-mp-1);
            end
            if ii<=n-mp*2
                residright(:,ii)=ptright*pp-x(idx+ii+mp-1);
            end
        end
        
        pleft=pt\x(idx);
        pright=pt\x(idx+n-mp);
        %
        paddedX=[ptleft*pleft;x;ptright*pright];
        %
        
        
    case {'minimumslope','minslope','ms'}
        idx=(1:mp)';
        
        residleft=nan(mp,n-mp*2);
        residright=nan(mp,n-mp*2);
        for ii=1:n-mp
            pp=mean(x(idx+ii-1));
            if ii>mp
                residleft(:,ii-mp)=pp-x(idx+ii-mp-1);
            end
            if ii<=n-mp*2
                residright(:,ii)=pp-x(idx+ii+mp-1);
            end
        end
        
        pleft=zeros(mp,1)+mean(x(idx));
        pright=zeros(mp,1)+mean(x(idx+n-mp));
        %
        paddedX=[pleft;x;pright];
        

end


if (mp*7>n)
    warning('The lag is quite big. That makes it difficult for the boostrapping procedure. Assuming symmetry of residual errors to get more bootstrap data...')
    
    residleft=[residleft flipud(residright)];
    residleft=[residleft -residleft];
    residright=flipud(residleft); 
end


% R=filter(tfilt,1,paddedX);
% R=R((1:n)+M*2-1);
R=filter2(tfilt,paddedX,'valid');

%What happens to whitenoise if it goes through this tfilt?
if noerrorbar
    sE_R=[];
else
    %pad sE_x
    %    sE_x=[zeros(mp,1)+epleft;sE_x;zeros(mp,1)+epright];
    
    %MONTE CARLO WAY
    mcreal=nan(length(x),2000);
    
    for ii=1:size(mcreal,2)
        noise=[residleft(:,ceil(rand*size(residleft,2)));sE_x();residright(:,ceil(rand*size(residright,2)))];
        mcreal(:,ii)=filter2(tfilt,paddedX+noise,'valid');
    end
    sE_R=std(mcreal')';
end

if (nargout==0)
    if ~isempty(sE_R)
        CI=prctile(mcreal',[5 95])';
        h=fill([t;flipud(t)],[CI(:,1);flipud(CI(:,2))],[.6 .6 .8]);
        set(h,'linestyle','none')
        hold on
    end
    hh=plot(t,x,'b',t,R,'k');
    set(hh(2),'linewidth',2);
    set(hh(1),'color',[.8 .8 1]);
    hold off
end


if nargout<1, clear R,end
if nargout<2, clear V,end
if nargout<3, clear sE_R,end
if nargout<4, clear trendidx,end
if nargout<5, clear tfilt,end

end


function ArgStruct=parseArgs(args,ArgStruct,varargin)
% Helper function for parsing varargin. 
%
%
% ArgStruct=parseArgs(varargin,ArgStruct[,FlagtypeParams[,Aliases]])
%
% * ArgStruct is the structure full of named arguments with default values.
% * Flagtype params is params that don't require a value. (the value will be set to 1 if it is present)
% * Aliases can be used to map one argument-name to several argstruct fields
%
%
% example usage: 
% --------------
% function parseargtest(varargin)
%
% %define the acceptable named arguments and assign default values
% Args=struct('Holdaxis',0, ...
%        'SpacingVertical',0.05,'SpacingHorizontal',0.05, ...
%        'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%        'MarginLeft',.1,'MarginRight',.1,'MarginTop',.1,'MarginBottom',.1, ...
%        'rows',[],'cols',[]); 
%
% %The capital letters define abrreviations.  
% %  Eg. parseargtest('spacingvertical',0) is equivalent to  parseargtest('sv',0) 
%
% Args=parseArgs(varargin,Args, ... % fill the arg-struct with values entered by the user
%           {'Holdaxis'}, ... %this argument has no value (flag-type)
%           {'Spacing' {'sh','sv'}; 'Padding' {'pl','pr','pt','pb'}; 'Margin' {'ml','mr','mt','mb'}});
%
% disp(Args)
%
%
%
%
% Aslak Grinsted 2004

% -------------------------------------------------------------------------
%   Copyright (C) 2002-2004, Aslak Grinsted
%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.

persistent matlabver

if isempty(matlabver)
    matlabver=ver('MATLAB');
    matlabver=str2double(matlabver.Version);
end

Aliases={};
FlagTypeParams='';

if (length(varargin)>0) 
    FlagTypeParams=lower(strvcat(varargin{1}));  %#ok
    if length(varargin)>1
        Aliases=varargin{2};
    end
end
 

%---------------Get "numeric" arguments
NumArgCount=1;
while (NumArgCount<=size(args,2))&&(~ischar(args{NumArgCount}))
    NumArgCount=NumArgCount+1;
end
NumArgCount=NumArgCount-1;
if (NumArgCount>0)
    ArgStruct.NumericArguments={args{1:NumArgCount}};
else
    ArgStruct.NumericArguments={};
end 


%--------------Make an accepted fieldname matrix (case insensitive)
Fnames=fieldnames(ArgStruct);
for i=1:length(Fnames)
    name=lower(Fnames{i,1});
    Fnames{i,2}=name; %col2=lower
    Fnames{i,3}=[name(Fnames{i,1}~=name) ' ']; %col3=abreviation letters (those that are uppercase in the ArgStruct) e.g. SpacingHoriz->sh
    %the space prevents strvcat from removing empty lines
    Fnames{i,4}=isempty(strmatch(Fnames{i,2},FlagTypeParams)); %Does this parameter have a value?
end
FnamesFull=strvcat(Fnames{:,2}); %#ok
FnamesAbbr=strvcat(Fnames{:,3}); %#ok

if length(Aliases)>0  
    for i=1:length(Aliases)
        name=lower(Aliases{i,1});
        FieldIdx=strmatch(name,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(name,FnamesFull); %&??????? exact or not? 
        end
        Aliases{i,2}=FieldIdx;
        Aliases{i,3}=[name(Aliases{i,1}~=name) ' ']; %the space prevents strvcat from removing empty lines
        Aliases{i,1}=name; %dont need the name in uppercase anymore for aliases
    end
    %Append aliases to the end of FnamesFull and FnamesAbbr
    FnamesFull=strvcat(FnamesFull,strvcat(Aliases{:,1})); %#ok
    FnamesAbbr=strvcat(FnamesAbbr,strvcat(Aliases{:,3})); %#ok
end

%--------------get parameters--------------------
l=NumArgCount+1; 
while (l<=length(args))
    a=args{l};
    if ischar(a)
        paramHasValue=1; % assume that the parameter has is of type 'param',value
        a=lower(a);
        FieldIdx=strmatch(a,FnamesAbbr,'exact'); %try abbreviations (must be exact)
        if isempty(FieldIdx) 
            FieldIdx=strmatch(a,FnamesFull); 
        end
        if (length(FieldIdx)>1) %shortest fieldname should win 
            [mx,mxi]=max(sum(FnamesFull(FieldIdx,:)==' ',2));%#ok
            FieldIdx=FieldIdx(mxi);
        end
        if FieldIdx>length(Fnames) %then it's an alias type.
            FieldIdx=Aliases{FieldIdx-length(Fnames),2}; 
        end
        
        if isempty(FieldIdx) 
            error(['Unknown named parameter: ' a])
        end
        for curField=FieldIdx' %if it is an alias it could be more than one.
            if (Fnames{curField,4})
                if (l+1>length(args))
                    error(['Expected a value for parameter: ' Fnames{curField,1}])
                end
                val=args{l+1};
            else %FLAG PARAMETER
                if (l<length(args)) %there might be a explicitly specified value for the flag
                    val=args{l+1};
                    if isnumeric(val)
                        if (numel(val)==1)
                            val=logical(val);
                        else
                            error(['Invalid value for flag-parameter: ' Fnames{curField,1}])
                        end
                    else
                        val=true;
                        paramHasValue=0; 
                    end
                else
                    val=true;
                    paramHasValue=0; 
                end
            end
            if matlabver>=6
                ArgStruct.(Fnames{curField,1})=val; %try the line below if you get an error here
            else
                ArgStruct=setfield(ArgStruct,Fnames{curField,1},val); %#ok <-works in old matlab versions
            end
        end
        l=l+1+paramHasValue; %if a wildcard matches more than one
    else
        error(['Expected a named parameter: ' num2str(a)])
    end
end

end


function [d,dt]=formatts(d)
%
% Usage: [d,dt]=formatts(d)
%
% Helper function for CWT,XWT,WTC
%
% Brings a timeseries into a shape so that it has two columns: [time, value].
%
%
% (C) Aslak Grinsted 2002-2014
% http://www.glaciology.net/wavelet-coherence

% -------------------------------------------------------------------------
%The MIT License (MIT)
%
%Copyright (c) 2014 Aslak Grinsted
%
%Permission is hereby granted, free of charge, to any person obtaining a copy
%of this software and associated documentation files (the "Software"), to deal
%in the Software without restriction, including without limitation the rights
%to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%copies of the Software, and to permit persons to whom the Software is
%furnished to do so, subject to the following conditions:
%
%The above copyright notice and this permission notice shall be included in
%all copies or substantial portions of the Software.
%
%THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
%THE SOFTWARE.
%---------------------------------------------------------------------------




if (ndims(d)>2)
    error('Input time series should be 2 dimensional.');
end
if (numel(d)==length(d))
    d=[(1:length(d))' d(:)];
end
if size(d,1)<size(d,2)
    d=d';
end
if (size(d,2)~=2)
    error('Time series must have 2 columns.')
end

if (d(2,1)-d(1,1)<0)
    d=flipud(d);
end

dt=diff(d(:,1));
if any(abs(dt-dt(1))>1e-1*dt(1))
    error('Time step must be constant.');
end
if (dt==0)
    error('Time step must be greater than zero.')
end

dt=dt(1);

end


function [E,V,A,R,C]=ssa(x,M,method) 
%  SSA - driver routine to perform Singular Spectrum Analysis
%  Syntax: [E,V,A,R,C]=ssa(x,M); [E,V,A,R,C]=ssa(x,M,'BK');
%
% Input:    x - time series
%           M - embedding dimension.
%      method - (optional) method of calculating the covariance matrix:
%                 'unbiased' (N-k weighted) (default)
%                 'biased'   (N-weighted or Yule-Walker)  
%                 'BK'       (Broomhead/King type estimate)
%
% Output: E - eigenfunction (T-EOF) matrix in standard form
%         V - vector containing variances (unnormalized eigenvalues)
%         A - Matrix of principal components
%         R - Matrix of reconstructed components
%         C - Covariance matrix
%
%  See Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
%
%  Written by Eric Breitenberger.     Version date 1/22/96
%  Please send comments and suggestions to eric@gi.alaska.edu   

if nargin==2, method='unbiased'; end

[E,V,C]=ssaeig(x,M,method);
[A]=pc(x,E);
[R]=rc(A,E);

end


  function [E,V,C]=ssaeig(x,M,method)
%  SSAEIG - starts an SSA of series 'x', for embedding dimension 'M'.
%  Syntax: [E,V,C]=ssaeig(x,M); [E,V,C]=ssaeig(x,M,'BK');  
%
% Input:    x - time series
%           M - embedding dimension.
%      method - (optional) method of calculating the covariance matrix:
%                 'unbiased' (N-k weighted) (default)
%                 'biased'   (N-weighted or Yule-Walker)  
%                 'BK'       (Broomhead/King type estimate)
%
% Output:   E - eigenfunction matrix in standard form
%               (columns are the eigenvectors, or T-EOFs)
%           V - vector containing variances (unnormalized eigenvalues)
%           C - covariance matrix
%
%  E and V are ordered from large to small.
%  See section 2 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
%
%  Written by Eric Breitenberger.    Version date 1/22/96
%  Please send comments and suggestions to eric@gi.alaska.edu   
%
if nargin==2, method='unbiased'; end
[N,col]=size(x);                
if col>N, x=x'; [N,col]=size(x); end   % change x to column vector if needed
if M>=N-M+1, error('Hey! Too big a lag!'), end
if col>=2, error('Hey! Vectors only!'), end

if ~strcmp(method,'BK')
  c=ac(x, M-1,method);   % calculate autocovariance estimates
  C=toeplitz(c);         % create Toeplitz matrix (trajectory matrix)
else 
  C=bk(x,M);           % Broomhead/King estimate
end

[E,L]=eig(C);          % calculate eigenvectors, values of C
[V,i]=sort(-diag(L));  % create sorted eigenvalue vector
V=-V';   
E=E(:,i);              % sort eigenvector matrix

  end
  

function  c=ac(x,k,method)
%    Syntax c=ac(x,k);  c=ac(x,k,'biased');  
% AC calculates the auto-covariance for series x out to lag k. The 
% result is output in c, which has k+1 elements. The first element
% is the covariance at lag zero; succeeding elements 2:k+1 are the
% covariances at lags 1 to k.
%
% Method can be: - 'biased'   (N-weighted or Yule-Walker)  
%                - 'unbiased' (N-k weighted - this is the default)
%                - 'BK'       (Broomhead/King type estimate)
%
% Note that the BK method is of limited use unless you are computing
% the full covariance matrix - see BK.M.                    
%
% Written by Eric Breitenberger.   Version date 1/11/96
% Please send comments and suggestions to eric@gi.alaska.edu   
 
[N,col]=size(x);
if col>N, tmp=col; col=N; N=tmp; end
if k>=N, error('Too big a lag!'), end
if col>1, error('Vectors only!'), end

x=x(:);
x=x-mean(x); % center the series.
c=zeros(1,k+1);

if nargin==2, method='unbiased'; end

if strcmp(method, 'unbiased')
  for i=1:k+1
   c(i)=(x(1:N-i+1)'*x(i:N))/(N-i+1);
  end

elseif strcmp(method,'biased')
  for i=1:k+1
    c(i)=x(1:N-i+1)'*x(i:N);
  end
  c=c./N;

elseif strcmp(method,'BK')
  for i=1:k+1
    c(i)=x(1:N-k)'*x(i:N-k+i-1);
  end
  c=c./(N-k);

else
  error('Improper specification of method.')
end

end


function [A]=pc(x, E)
% PC - calculates principal components
%    Syntax: [A]=pc(x, E); 
%  PC calculates the principal components of the series x
%  from the eigenfunction matrix E (output from SSAEIG).
%  Returns:      A - principal components matrix (N-M+1 x M)
%  See section 2.4 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
%
%  Written by Eric Breitenberger.    Version date 5/20/95
%  Please send comments and suggestions to eric@gi.alaska.edu   
%
[N,col]=size(x);
if min(N,col)>1, error('x must be a vector.'), end
if col>1, x=x'; N=col; end     % convert x to column if necessary.
x=x-mean(x);

[M,c]=size(E);                
if M~=c, error('E is improperly dimensioned'), end
A=zeros(N-M+1,M);
% This could be rewritten using 'filter', like MPC .
for i=1:N-M+1;                 
  w=x(i:i+M-1);          
  A(i,:)=w'*E;
end

end


function [R]=rc(A,E)
% Syntax: [R]=rc(A,E);
% This function calculates the 'reconstructed components' using the 
% eigenvectors (E, from ssaeig.m) and principal components (A, from pc.m).
% R is N x M, where M is the embedding dimension used in ssaeig.m.
%
% See section 2.5 of Vautard, Yiou, and Ghil, Physica D 58, 95-126, 1992.
%
%  Written by Eric Breitenberger.   Version date 5/18/95
%  Please send comments and suggestions to eric@gi.alaska.edu   

[M,c]=size(E);
[ra, ca]=size(A);
if M~=c, error('E is improperly dimensioned.'),end
if ca~=M, error('A is improperly dimensioned.'),end
N=ra+M-1;  % Assumes A has N-M+1 rows.

R=zeros(N,M);
Z=zeros(M-1,M);
A=[A' Z'];
A=A';

% Calculate RCs
for k=1:M
  R(:,k)=filter(E(:,k),M,A(:,k));
end

% Adjust first M-1 rows and last M-1 rows
for i=1:M-1
  R(i,:)=R(i,:)*(M/i);
  R(N-i+1,:)=R(N-i+1,:)*(M/i);
end

end
  


