function [cOptLLR, cOptLLRCons, output] = c_opt_llr(x, y, level, show, h_opt)
% FUNCTION [COPTLLR, COPTLLRCONS, OUTPUT] = C_OPT_LLR(X, Y, LEVEL, SHOW,
% H_OPT)
% Returns optimal thresholds for an RDD
% =======================================================================
% Inputs
% =======================================================================
% x                       = running variable (discontinuity at x=0)
% y                       = outcome variable
% level                   = level of significance for conservative
%                           threshold, default  .05
% show                    = display data and fit regressions, default F
% h_opt                   = bandwidth for data used in the regression,
%                           default value is optimal IK bandwidth with 
%                           uniform kernel
% =======================================================================
% Output
% =======================================================================
% cOptLLR                = optimal threshold 
% cOptLLRCons            = conservative optimal threshold
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
if nargin==4
[~,~,h_opt]=rd_optbandwidth_uni(y, x, x>0,[],0,false,0);
end

output.h_opt= h_opt;
indsIK = x>-h_opt & x<h_opt; % keep only values within bandwidht
xIK = x(indsIK); yIK = y(indsIK);
xIKL = xIK(xIK<0); xIKR = xIK(xIK>0); % split left and right data set
yIKL = yIK(xIK<0); yIKR = yIK(xIK>0); 


[llrMdl, llrMdlStats] = polyfit(xIKL, yIKL, 1);
[llrMdr, llrMdrStats] = polyfit(xIKR, yIKR, 1);
output.modelL= llrMdl; output.modelR = llrMdr;

% LATE estimate
[~, output.se, output.tau] = tau_low(0,  llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
% standard errors and p-value
z = abs(output.tau/output.se);
output.pval =  exp(-0.717*z - 0.416*(z^2));
% confIntL = conf_int_llr(-h_opt, h_opt, xIKL, yIKL, llrMdl, level);
% confIntR =  conf_int_llr(-h_opt, h_opt, xIKR, yIKR, llrMdr, level);

% optimal threshold
tmp =llrMdr-llrMdl;
cOptLLR =- tmp(2)/tmp(1);

if cOptLLR>h_opt;cOptLLR=h_opt;end % keep inside bw
if cOptLLR<-h_opt;cOptLLR=-h_opt;end

% conservative threshold
fun = @(x) tau_low(x, llrMdl, llrMdr, xIKL, xIKR, yIKL, yIKR, level);
[cOptLLRCons, ~, exitflag]  = fzero(fun, 0);
if exitflag~=1
    cOptLLRCons = nan;
elseif cOptLLRCons>h_opt
    cOptLLRCons=h_opt;% keep inside bw
elseif cOptLLRCons<-h_opt
    cOptLLRCons=-h_opt;
end



if show    % show data, bw and regression lines
    confIntL = conf_int_llr(linspace(-h_opt, h_opt), xIKL, yIKL, llrMdl, level);
    confIntR =  conf_int_llr(linspace(-h_opt, h_opt), xIKR, yIKR, llrMdr, level);
    cm = jet;
    figure(11); ah = gca; 
    legH =[]; % legend handles
    % plot data

    l=scatter(x, y,  'Marker', '.', 'DisplayName', 'data');legH = [legH l];
    % hf = fill_between(xmean,ymean-ystd,ymean+ystd);
    % hf.FaceColor = [.9 .9 .9];
    % 
    hold on%%

    l =line([0, 0], ah.YLim, 'linewidth', 1, 'color', 'k', 'DisplayName','c');legH = [legH l];
    % plot KI optimal bandwidth
    l =line([h_opt, h_opt], ah.YLim,'linestyle',  '--', 'linewidth', 1, 'color', 'k', 'DisplayName','IK');legH = [legH l];
    line([-h_opt, -h_opt], ah.YLim, 'linestyle', '--', 'linewidth', 1, 'color', 'k');
    % plot LLR extrapolation
    l=plot([-h_opt, 0], polyval(llrMdl,[-h_opt, 0]), 'linewidth', 1.5, 'color' , cm(1, :), 'DisplayName','LLR');legH = [legH l];
    plot([0, h_opt], polyval(llrMdr,[0, 0+h_opt]),'linewidth',1.5, 'color', cm(1, :));
    plot([0-h_opt, 0], polyval(llrMdr,[0-h_opt, 0]), 'linewidth', 1.5, 'color' , cm(1, :), 'linestyle', '-.');

    hf = fill_between(linspace(-h_opt, h_opt), confIntL(:, 1), confIntL(:, 2));hf.FaceColor = [.9 .9 .9];
    hf = fill_between(linspace(-h_opt, h_opt), confIntR(:, 1), confIntR(:, 2));hf.FaceColor = [.9 .9 .9];
end

end

function [tau_low, se, tau] =  tau_low(x_new,  modelL, modelR, x_l, x_r, y_l, y_r, level)
% estimates late, standard erros and conservative late at x=x_new
% (c) sof.triantafillou@gmail.com
tau = polyval(modelR, x_new)-polyval(modelL, x_new);

var_l = var_pred(x_new, x_l,  y_l, polyval(modelL, x_l));
var_r = var_pred(x_new, x_r,  y_r, polyval(modelR, x_r));

se = sqrt(var_l+var_r)'; 
quant = -norminv(level);
tau_low = tau -se*quant;
end

function var_pred = var_pred(x_new, x, y, y_hat)
% conditional variance var(x_new|y)
% (c) sof.triantafillou@gmail.com
n = length(y); % sample size
x_new= reshape(x_new, length(x_new),1); %make column vector;
eps = y_hat-y; 
s2_yx = sum(eps.^2)/(n-2);

x_new_= [x_new ones(length(x_new),1)];
x_ = [x ones(length(x), 1)];

var_pred = nan(size(x_new, 1), 1);
%\hat V_p=s^2 \cdot \mathbf{x_0} \cdot(\mathbf{X'X})^{-1}\mathbf{x'_0},'*x
for i=1:size(x_new, 1)
    var_pred(i) = s2_yx*(x_new_(i, :)*inv(x_'*x_)*x_new_(i, :)');
end

end

