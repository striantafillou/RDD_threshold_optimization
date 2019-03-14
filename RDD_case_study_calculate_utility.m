clear;
clc;
%control panel
load dftips.mat
level = 0.05;

dftips = dftips(dftips.pkp_day.Month==1, :); % keep january data only
x = dftips.fare-15;y = dftips.tip_frac;% recenter around theshold
%x = x(y<200); y = y(y<200);

%% LLR with uniform kernel in IK bandwidth
[est,se,h_opt]=rd_optbandwidth_uni(y, x, x>0,[],0,false,0);

indsIK = x>-h_opt & x<h_opt; % keep only values within bandwidth
xIK = x(indsIK); yIK = y(indsIK);
xIKL = xIK(xIK<0); xIKR = xIK(xIK>0); % split left and right data set
yIKL = yIK(xIK<0); yIKR = yIK(xIK>0); 


[llrMdl, llrMdlStats] = polyfit(xIKL, yIKL, 1);
[llrMdr, llrMdrStats] = polyfit(xIKR, yIKR, 1);


confIntL = conf_int_llr(-h_opt, h_opt, xIKL, yIKL, llrMdl, level);
confIntR =  conf_int_llr(-h_opt, h_opt, xIKR, yIKR, llrMdr, level);


tmp =llrMdr-llrMdl;
cLLR =- tmp(2)/tmp(1);


fun = @(x) tau_low(x, llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
cLLRCons = fzero(fun, 0);
%% calculate as discrete sum

%clc;clear;load casestudy.mat
% calculate mean value of f_l, f_r in the interval (0, cLLR)
[xmean, ymean, ystd, nPoints]  = meanbin(x, y);% desnsity of observations is assumed uniform in [cLLR, 0] with density

xLs = xmean(xmean<0 & xmean>cLLR)+15;
nPointsL = nPoints(xmean<0 & xmean>cLLR);
xDensL = nPointsL/sum(nPointsL);
%xDensL = xDens(xmean<0& xmean>cLLR);

% % recenter around 15
% llrMdl_new = [llrMdl(1) llrMdl(2)- 15*llrMdl(1)];
% llrMdr_new = [llrMdr(1) llrMdr(2)- 15*llrMdr(1)];
% yLs = polyval(llrMdl_new, xLs);
% yLs_extrap = polyval(llrMdr_new, xLs);

% calculate tip percentages.
yLs = polyval(llrMdl, xLs-15);
yLs_extrap = polyval(llrMdr, xLs-15);

percL = sum(yLs.*xDensL);
percL_extrap = sum(yLs_extrap.*xDensL);

perc_gain = percL_extrap- percL


tipL= sum(yLs.*xDensL.*xLs)/100;
tipL_extrap = sum(yLs_extrap.*xDensL.*xLs)/100;

tip_gain = tipL_extrap-tipL
money_gain = tip_gain*sum(nPointsL)

%% for conservative estimate

consInds = xLs>(cLLRCons+15);
xLsCons = xLs(consInds);
nPointsLCons = nPointsL(consInds);
xDensLCons =  nPointsLCons/sum(nPointsLCons);
%xDensL = xDens(xmean<0& xmean>cLLR);

% calculate tip percentages.
yLsCons = polyval(llrMdl, xLsCons-15);
yLs_extrapCons = polyval(llrMdr, xLsCons-15);

percLCons = sum(yLsCons.*xDensLCons);
percL_extrapCons = sum(yLs_extrapCons.*xDensLCons);

perc_gainCons = percL_extrapCons- percLCons


tipLCons= sum(yLsCons.*xDensLCons.*xLsCons)/100;
tipL_extrapCons = sum(yLs_extrapCons.*xDensLCons.*xLsCons)/100;

tip_gainCons = tipL_extrapCons-tipLCons
money_gainCons = tip_gainCons*sum(nPointsLCons)





%%
clc;clear;load casestudy.mat
% calculate mean value of f_l, f_r in the interval (0, cLLR)
syms x
% desnsity of observations is assumed uniform in [cLLR, 0] with density
% 1/-cLLR
f_l =(llrMdl(1)*(x-15)+llrMdl(2))/100; % divide by 100 to get percentage, recenter around 15.
f_r =(llrMdr(1)*(x-15)+llrMdr(2))/100;
a1 =int(f_r*x*((-cLLR)^-1), x,  cLLR+15, 0+15);% integral f_x*x*p(x).
b1 =int(f_l*x*((-cLLR)^-1), x,  cLLR+15, 0+15);
% take difference
t1 =double(a1-b1);
%
a2 =int(f_r*x*((-cLLR)^-1), x,  cLLRCons+15, 0+15);
b2 =int(f_l*x*((-cLLR)^-1), x,  cLLRCons+15, 0+15);
t2 =double(a2-b2)
% take difference


%% calculate triangle area 
Ax = cLLR; Ay = polyval(llrMdr, cLLR);
Bx = 0; By = polyval(llrMdr, 0);
Cx = 0; Cy = polyval(llrMdl, 0);

area = abs(0.5*	(Ax*(By-Cy)+Bx*(Cy-Ay)+ Cx*(Ay-By)));
%% calculate triangle area for conservative 
Ax = cLLR; Ay = polyval(llrMdr, cLLR);
Bx = cLLRCons; By = polyval(llrMdr, cLLRCons);
Cx = cLLRCons; Cy = polyval(llrMdl,cLLRCons);
area_small=abs(0.5*	(Ax*(By-Cy)+Bx*(Cy-Ay)+ Cx*(Ay-By)));
area2= area-area_small