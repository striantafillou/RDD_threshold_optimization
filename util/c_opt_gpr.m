function [cOptGPR, cOptGPRCons, output] = c_opt_gpr(x, y, level, show)
% FUNCTION [COPTGPR, COPTGPRCONS, OUTPUT] = C_OPT_LLR(X, Y, LEVEL, SHOW)
% Returns optimal thresholds for an RDD
% =======================================================================
% Input
% =======================================================================
% x                       = running variable (discontinuity at x=0)
% y                       = outcome variable
% level                   = level of significance for conservative
%                           threshold, default  .05
% show                    = display data and fit regressions, default F
% =======================================================================
% Output
% =======================================================================
% cOptGPR                = optimal threshold 
% cOptGPRCons            = conservative optimal threshold
% output
%   .modelL, modelR      = left and right regression models
%   .tau                 = LATE at x=0
%   .se                  = standard error for LATE
%   .h_opt               = bandwidth used
%   .pval                = pvalue for null hypothesis of LATE=0
% (c) sof.triantafillou@gmail.com

if nargin==2
    level=0.05;show=false;
elseif nargin==3
    show=false;
end

xL = x(x<0); xR = x(x>0); % split left and right data set
yL = y(x<0); yR = y(x>0); 

% fit GP regressions

gprMdlr = fitrgp(xR, yR, 'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus', 'ShowPlots', false', 'verbose', 0));

gprMdll = fitrgp(xL, yL, 'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus', 'ShowPlots', false', 'verbose', 0));
output.modelL= gprMdll; output.modelR = gprMdlr;

% LATE estimate
[y_r, se_r] = predict(gprMdlr, 0);
var_r = se_r^2;
[y_l, se_l] = predict(gprMdll, 0);
var_l = se_l^2;
output.tau = y_r-y_l;
% standard error and p-value
output.se = sqrt(var_l+var_r);
z = abs(output.tau/output.se);
output.pval =  exp(-0.717*z - 0.416*(z^2));

% optimal threshold
fun1 = @(x) (predict(gprMdlr, x)-predict(gprMdll, x));
[cOptGPR, ~, exitflag]  = fzero(fun1, 0);
if exitflag~=1
    cOptGPR = nan;
elseif cOptGPR>max(x)
    cOptGPR=max(x);% keep inside bw
elseif cOptGPR<min(x)
    cOptGPR=min(x);
end

% conservative threshold
fun = @(x) tau_low_gpr(x, gprMdll, gprMdlr, level);
[cOptGPRCons, ~, exitflag]  = fzero(fun, 0);
if exitflag~=1
    cOptGPRCons = nan;
elseif cOptGPRCons>max(x)
    cOptGPRCons=max(x);% keep inside bw
elseif cOptGPRCons<min(x)
    cOptGPRCons=min(x);
end


if show
    cm = brewermap(10, 'set2');
    figure(22); ah = gca; 
    legH =[]; % legend handles
    % plot data

    l=scatter(x, y,  'Marker', '.', 'DisplayName', 'data');legH = [legH l];
    % hf = fill_between(xmean,ymean-ystd,ymean+ystd);
    % hf.FaceColor = [.9 .9 .9];
    % 
    hold on%%
 
    l =line([0, 0], ah.YLim, 'linewidth', 1, 'color', 'k', 'DisplayName','c');legH = [legH l];
    % plot  extrapolation
    xPredL = linspace(-1, 0)';
    xPredR = linspace(0, 1)';
    [yPredL, ~, yPredLCI] = predict(gprMdll, xPredL);
    [yPredR, ~, yPredRCI] = predict(gprMdlr, xPredR);
    [yPredRE, ~, yPredRECI] = predict(gprMdlr, xPredL);
    l=plot(xPredL, yPredL, 'linewidth', 1.5, 'color' , cm(1, :), 'DisplayName','GPR');legH = [legH l];
    plot(xPredR, yPredR,'linewidth',1.5, 'color', cm(1, :));
    plot(xPredL, yPredRE, 'linewidth', 1.5, 'color' , cm(1, :), 'linestyle', '-.');
    

    hf = fill_between(xPredL, yPredLCI(:, 1), yPredLCI(:, 2));hf.FaceColor = [.9 .9 .9];
    hf = fill_between(xPredR, yPredRCI(:, 1), yPredRCI(:, 2));hf.FaceColor = [.9 .9 .9];
    hf = fill_between(xPredL, yPredRECI(:, 1), yPredRECI(:, 2));hf.FaceColor = [.9 .9 .9];

end

end

function [tau_low] =  tau_low_gpr(x_new,  modelL, modelR,  level)
% estimates late, standard erros and conservative late at x=x_new
% (c) sof.triantafillou@gmail.com
[y_r, se_r] = predict(modelR, x_new);
var_r = se_r^2;
[y_l, se_l] = predict(modelL, x_new);
var_l = se_l^2;
tau = y_r-y_l;



se = sqrt(var_l+var_r);
quant = -norminv(level);
tau_low = tau -se*quant;
end

