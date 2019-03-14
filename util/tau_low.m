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


