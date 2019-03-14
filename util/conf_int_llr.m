function confInt = conf_int_llr(x_new, x, y, model, level)
% calculates confidence intervals for y=ax+b at x_new
x_new= reshape(x_new, length(x_new),1);
var = var_pred(x_new, x, y, polyval(model, x)); % variance 
% confidence intervals
quant = -norminv(level/2);
y_new = polyval(model, x_new);
confInt = [y_new- quant*sqrt(var) y_new + quant*sqrt(var)];
end