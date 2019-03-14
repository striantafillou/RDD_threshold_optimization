function [est,se,h_opt]=rd_optbandwidth_uni(y,x,w,z,threshold,output,addcov)

% calculates optimal IK bandwidth for regression discontinuity settings

% INPUTS
% y is the outcome
% x is the scalar forcing variable
% w is the binary treatment indicator (equal to an indicator that x exceeds threshold in a sharp 
%    rd design.
% z is a vector of additional covariates (always need as argument, even if not used)
% threshold is the value of the threshold for the forcing variable
% if output=1, print all intermediate steps, if output=0 do not
% if addcov=1, the variables in z will be used as additional covariates
%    if addcov=0 they will be ignored.

% OUTPUT	
% h_opt is bandwidth
% est is point estimate
% se is standard error

c=threshold;
n=length(y);    % number of observations
x=x-c;          % forcing variable relative to threshold
mx=mean(x);
sx=sqrt((x-mx)'*(x-mx)/(n-1));  % sample standard deviation of forcing variable

%pause

h_silverman=1.84*sx*(n^(-1/5));     % initial bandwidth
                                   % based on silverman's rule

%fprintf('%.3f\n', h_silverman)
n_min=sum(x<0);           % number of observations to left of threshold
n_plus=sum(x>=0);         % number of observations to the right of threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% optimal bandwidth                                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ===============================================================================
% step one

i_plus_1=(x>=0)&(x<=h_silverman);   % positive observatins within bandwidth
i_min_1=(x<0)&(x>=-h_silverman);    % negative observations within bandwidth

n_plus_1=sum(i_plus_1);    % number of positive observations within bandwidth
n_min_1=sum(i_min_1);      % number of negative observations within bandwidth

y_ave_plus_1=y'*i_plus_1/n_plus_1;   % average outcome on positive side within bandwidth
y_ave_min_1=y'*i_min_1/n_min_1;      % average outcome on negative side within bandwidth

dy_plus_1=y(i_plus_1,1)-y_ave_plus_1; % residuals on positive side within bandwidth
dy_min_1=y(i_min_1,1)-y_ave_min_1;    % residuals on negative side within bandwidth

s_min_1=sqrt(dy_min_1'*dy_min_1/(n_min_1-1));     % standard dev for obs to left of threshold within bandwidth
s_plus_1=sqrt(dy_plus_1'*dy_plus_1/(n_plus_1-1)); % standard dev for obs to right of threshold within bandwidth

% estimated variance at threshold
sigmas=(dy_plus_1'*dy_plus_1+dy_min_1'*dy_min_1)/(n_plus_1+n_min_1);
% estimated standard deviation at threshold
sigma=sqrt(sigmas);            

% estimated density at threshold
fc=(n_plus_1+n_min_1)/(2*n*h_silverman);  % estimate of marginal density at threshold





% ===============================================================================
% step two

median_plus=median(x(x>=0,1));
median_min=median(x(x<0,1));
middle=(x>=median_min)&(x<=median_plus);
x_middle=x(middle,1);
y_middle=y(middle,1);
tt=[ones(sum(middle),1),x_middle>=0,x_middle,x_middle.*x_middle,x_middle.*x_middle.*x_middle];
%tt=[[t_min,zeros(n_min_2,1)];[t_plus,ones(n_plus_2,1)]];
gamma=inv(tt'*tt)*(tt'*y_middle);
yhat=tt*gamma;
third_der=6*gamma(5,1);           % estimate of third derivative

% optimal bandwidth for estimating second derivatives
h_plus_2=3.56*(sigmas/(fc*max(third_der*third_der,0.01)))^(1/7)*n_plus^(-1/7);
h_min_2=3.56*(sigmas/(fc*max(third_der*third_der,0.01)))^(1/7)*n_min^(-1/7);

n_plus_2=sum(middle&(x>=0));
n_min_2=sum(middle&(x<0));

i_min_3=(x<0)&(x>=-h_min_2);     % negative values within bandwidth
if sum(i_min_3)<=10,
   i_min_3=(x<0);
   end
n_min_3=sum(i_min_3);
y_min_3=y(i_min_3,1);            % y-values for negative x within bandwidth
x_min_3=x(i_min_3,1);            % x-values for negative x within bandwidth 
t_min=[ones(sum(i_min_3),1),x_min_3,x_min_3.*x_min_3];
% regression from the left
lambda_min=inv(t_min'*t_min)*(t_min'*y_min_3);
%' second derivative from left'
second_der_min=2*lambda_min(3,1);        % estimate of second derivative from left

i_plus_3=(x>=0)&(x<=h_plus_2);
n_plus_3=sum(i_plus_3);
if sum(i_plus_3)<=10,
   i_plus_3=(x>=0);
   end
y_plus_3=y(i_plus_3,1);
x_plus_3=x(i_plus_3,1);
t_plus=[ones(sum(i_plus_3),1),x_plus_3,x_plus_3.*x_plus_3];
%'lambda estimate from right'
lambda_plus=inv(t_plus'*t_plus)*(t_plus'*y_plus_3);
%' second derivative from right'
second_der_plus=2*lambda_plus(3,1);        % estimate of second derivative from right



% ===============================================================================
%'begin step 3'

%'r_plus regularization term from right'
r_plus=720*sigmas/(n_plus_3*(h_plus_2^4));
%'r_min regularization term from left'
r_min=720*sigmas/(n_min_3*(h_min_2^4));

CK=5.4;
%'hat h opt'
h_opt=CK*((2*sigmas/(fc*((second_der_plus-second_der_min)^2+r_plus+r_min)))^(1/5))*(n^(-1/5));



% ===============================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% estimation and inference                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weight=((abs(x)<=h_opt));   % weight for local linear regression
summ_weight=[mean(weight),std(weight)];
if addcov==1,
meanz=(z'*weight/sum(weight))';
z=z-ones(length(z),1)*meanz;                  % z in deviations from weighted mean
end
x_plus=x(x>=0,1);
y_plus=y(x>=0,1);
x_min=x(x<0,1);
y_min=y(x<0,1);
w_plus=w(x>=0,1);
w_min=w(x<0,1);
if addcov==1,
   z_plus=z(x>=0,:);
   z_min=z(x<0,:);
   xz_min=[x_min,z_min];
   xz_plus=[x_plus,z_plus];
   else
   xz_min=[x_min];
   xz_plus=[x_plus];
   end
[NN,KXZ]=size(xz_min);




n_h=sum(abs(x_plus)<h_opt);  % number of observations within bandwidth.
if n_h<5,      % too few observations close by,
               % and so we just average the closest five observations
   sort_x=sort(x_plus);
   limit=sort_x(min(5,length(x_plus)),1);
   mu_plus=mean(y_plus(x_plus<=limit,1));
   var_plus=(std(y_plus(x_plus<=limit,1))^2)/sum(x_plus<=limit,1);
   
   else,       % with a sufficient number of observations close by
               % we do local linear regression
   weight=(x_plus<=h_opt);%.*(1-x_plus/h_opt);     % weight for local linear regression
   sweight=sqrt(double(weight((x_plus<=h_opt),1)));      % square root of weights 
   xx_plus=[sweight,xz_plus(abs(x_plus)<=h_opt,:).*(sweight*ones(1,KXZ))];
   yy_plus=y_plus(x_plus<=h_opt).*sweight;
   beta_y_plus=inv(xx_plus'*xx_plus)*(xx_plus'*yy_plus);       % weighted regression
   mu_y_plus=beta_y_plus(1,1);                    % prediction at zero
   eps_y_plus=y_plus-[ones(length(y_plus),1),xz_plus]*beta_y_plus;
   xxx_y_plus=[eps_y_plus(x_plus<=h_opt,1).*weight(x_plus<=h_opt,1),xz_plus(x_plus<=h_opt,:).*((eps_y_plus(x_plus<=h_opt,1).*weight(x_plus<=h_opt,1))*ones(1,KXZ))];
   sig=eps_y_plus'*(weight.*eps_y_plus)/(length(y_plus)-2);
   var=inv(xx_plus'*xx_plus)*(xxx_y_plus'*xxx_y_plus)*inv(xx_plus'*xx_plus);
   delta_y_plus=(xxx_y_plus'*xxx_y_plus);
   gamma_plus=(xx_plus'*xx_plus);
   v_plus=var;
   var_plus=var(1,1);

   ww_plus=w_plus(x_plus<=h_opt).*sweight;
   beta_w_plus=inv(xx_plus'*xx_plus)*(xx_plus'*ww_plus);       % weighted regression
   mu_w_plus=beta_w_plus(1,1);                    % prediction at zero
   eps_w_plus=w_plus-[ones(length(w_plus),1),xz_plus]*beta_w_plus;
   xxx_w_plus=[eps_w_plus(x_plus<=h_opt,1).*weight(x_plus<=h_opt,1),xz_plus(x_plus<=h_opt,:).*((eps_w_plus(x_plus<=h_opt,1).*weight(x_plus<=h_opt,1))*ones(1,KXZ))];
   sig=eps_w_plus'*(weight.*eps_w_plus)/(length(w_plus)-2);
   var=inv(xx_plus'*xx_plus)*(xxx_y_plus'*xxx_y_plus)*inv(xx_plus'*xx_plus);
   delta_w_plus=(xxx_w_plus'*xxx_w_plus);
   delta_yw_plus=(xxx_y_plus'*xxx_w_plus);
   delta_yw_plus-delta_yw_plus';
   end

n_h=sum(abs(x_min)<h_opt);  % number of observations within bandwidth.
if n_h<5,      % too few observations close by,
               % and so we just average the closest five observations
   sort_x=sort(x_min);
   limit=sort_x(min(5,length(x_min)),1);
   mu_min=mean(y_min(x_min<=limit,1));
   var_min=(std(y_min(x_min<=limit,1))^2)/sum(x_min<=limit,1);
   
   else,       % with a sufficient number of observations close by
               % we do local linear regression
   weight=((abs(x_min)<=h_opt));%.*(1-abs(x_min)/h_opt);   % weight for local linear regression
   sweight=sqrt(double(weight(((abs(x_min)<=h_opt)),1)));      % square root of weights 
   xx_min=[sweight,xz_min(abs(x_min)<=h_opt,:).*(sweight*ones(1,KXZ))];
   yy_min=y_min((abs(x_min)<=h_opt)).*sweight;
   beta_y_min=inv(xx_min'*xx_min)*(xx_min'*yy_min);       % weighted regression
   mu_y_min=beta_y_min(1,1);                    % prediction at zero
   eps_y_min=y_min-[ones(length(y_min),1),xz_min]*beta_y_min;
   xxx_y_min=[eps_y_min(abs(x_min)<=h_opt,1).*weight(abs(x_min)<=h_opt,1),xz_min(abs(x_min)<=h_opt,:).*((eps_y_min(abs(x_min)<=h_opt,1).*weight(abs(x_min)<=h_opt,1))*ones(1,KXZ))];
   sig=eps_y_min'*(weight.*eps_y_min)/(length(y_min)-2);
   var=inv(xx_min'*xx_min)*(xxx_y_min'*xxx_y_min)*inv(xx_min'*xx_min);
   delta_y_min=(xxx_y_min'*xxx_y_min);
   gamma_min=(xx_min'*xx_min);
   v_min=var;
   var_min=var(1,1);
   ww_min=w_min((abs(x_min)<=h_opt)).*sweight;
   beta_w_min=inv(xx_min'*xx_min)*(xx_min'*ww_min);       % weighted regression
   mu_w_min=beta_w_min(1,1);                    % prediction at zero
   eps_w_min=w_min-[ones(length(w_min),1),xz_min]*beta_w_min;
   xxx_w_min=[eps_w_min(abs(x_min)<=h_opt,1).*weight(abs(x_min)<=h_opt,1),xz_min(abs(x_min)<=h_opt,:).*((eps_w_min(abs(x_min)<=h_opt,1).*weight(abs(x_min)<=h_opt,1))*ones(1,KXZ))];
   sig=eps_w_min'*(weight.*eps_w_min)/(length(w_min)-2);
   var=inv(xx_min'*xx_min)*(xxx_y_min'*xxx_y_min)*inv(xx_min'*xx_min);
   delta_w_min=(xxx_w_min'*xxx_w_min);
   delta_yw_min=(xxx_y_min'*xxx_w_min);
   end

est=(mu_y_plus-mu_y_min)/(mu_w_plus-mu_w_min);


% ===============================================================================

zzero=zeros(KXZ+1,KXZ+1);
delta=[[delta_y_min,delta_yw_min,zzero,zzero];[delta_yw_min',delta_w_min,zzero,zzero];[zzero,zzero,delta_y_plus,delta_yw_plus];[zzero,zzero,delta_yw_plus',delta_w_plus]];
gamma=[[gamma_min,zzero,zzero,zzero];[zzero,gamma_min,zzero,zzero];[zzero,zzero,gamma_plus,zzero];[zzero,zzero,zzero,gamma_plus]];

var=inv(gamma)'*delta*inv(gamma);
g=zeros(4*(KXZ+1),1);
g(0*(KXZ+1)+1,1)=-1/(mu_w_plus-mu_w_min);
g(1*(KXZ+1)+1,1)=(mu_y_plus-mu_y_min)/((mu_w_plus-mu_w_min)^2);
g(2*(KXZ+1)+1,1)=1/(mu_w_plus-mu_w_min);
g(3*(KXZ+1)+1,1)=-(mu_y_plus-mu_y_min)/((mu_w_plus-mu_w_min)^2);



vvv=g'*var*g;
svvv=sqrt(vvv);

if output==1,
   '==================================================================' %#ok<NOPRT>
%    'sample size, number of obs with x>=c, number of obs with x<c'
%    [n,n_plus,n_min]
%    'mean and standard deviation of forcing variable' 
%    [mean(x)+c,sx]
%    'silverman bandwidth'
%    h_silverman
%    'estimated standard deviation of outcome at threshold'
%    sigma
%    'estimated density at threshold'
%     fc   
%    'median of forcing variable to the left and right of the threshold'
%    [median_min,median_plus]+c
%    'coefficient in cubic regression'
%    third_der/6
%    'estimate of third derivative'
%    third_der
%    'bandwidth two'
%    [h_min_2,h_plus_2]
%    [c-h_min_2,c]
%    [c,c+h_plus_2]
%    'sample sizes to left and right within median'
%    [n_min_2,n_plus_2]
%    'second derivatives from left and right'
%    [second_der_min,second_der_plus]
%    'regularization terms'
%    [r_min,r_plus]
%    h_opt
%    summ_weight
%   meanz
   betas=[beta_y_min,beta_y_plus,beta_w_min,beta_w_plus]
%    delta_y_min
%    delta_w_min
%    delta_yw_min
%    gamma_min
%    v_min
%    se_min=sqrt(diag(v_min))
%    delta_y_plus
%    delta_w_plus
%    delta_yw_plus
%    gamma_plus
%    v_plus
%    se_plus=sqrt(diag(v_plus))
%    v_plus+v_min
%    delta
%    gamma
%    var
%    g
%    var_y_min=inv(gamma_min)*delta_y_min*inv(gamma_min)
%    se_y_min=sqrt(diag(var_y_min))
%    var_y_plus=inv(gamma_plus)*delta_y_plus*inv(gamma_plus)
%    se_y_plus=sqrt(diag(var_y_plus))
%    var_w_min=inv(gamma_min)*delta_w_min*inv(gamma_min)
%    se_w_min=sqrt(diag(var_w_min))
%    var_w_plus=inv(gamma_plus)*delta_w_plus*inv(gamma_plus)
%    se_w_plus=sqrt(diag(var_w_plus))
%   '==================================================================' 
   end


var_tau=g'*var*g;
tau=est;
se_tau=sqrt(var_tau);
se=se_tau;