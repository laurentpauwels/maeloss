%Choose Data file: Inflation rate, Growth rate, Unemployment rate
load('InflationData')
%load('GrowthData')
%load('UnempData')

y = actual;
Z = [fcast_95,fcast_94,fcast_37,fcast_89];
[T,p] = size(Z);
T2 = floor(T/2):T;%T2=floor(T/2);%
a_MSE = zeros(numel(T2),p);
a_MAE = zeros(numel(T2),p);
options = optimset('Display','notify-detailed','LargeScale','off',...
    'MaxFunEvals',5000,'MaxIter',1000);%'TolFun',1e-20,'TolX',1e-20);
for i = 1:numel(T2)
    t = T2(1,i);
    %MSE weights
    errors = (y(1:t,:).*ones(t,p))-Z(1:t,:);
    
    %Omegahat = (errors'*errors)/t;
    Omegahat = cov(errors);
    a_MSE(i,:) = Omegahat\ones(p,1)/(ones(p,1)'/Omegahat*ones(p,1));
    
    %MAE weights
    a_init = Omegahat\ones(p,1)/(ones(p,1)'/Omegahat*ones(p,1));
    f_mae = @(a_mae)maeloss(a_mae,errors);
    a_mae_raw = fminunc(f_mae, a_init(1:p-1,1),options);
    a_MAE(i,:) = [a_mae_raw;1-sum(a_mae_raw)];   
end

% Check Kurtosis and Skewness in individual forecast errors
ferrors = y(1:end,:) - Z(1:end,:);
errors_kurt = kurtosis(ferrors);
errors_skew = skewness(ferrors);

%% Table 1
%JB test for normality in individual forecast errors
jbpval = zeros(size(ferrors,2),1);
for i = 1:size(ferrors,2)
jb = (T/6)*((skewness(ferrors(:,i)).^2)+((kurtosis(ferrors(:,i)) - 3).^2) / 4);
jbpval(i,1) = 1-chi2cdf(jb,2);%Jarque-Berra test
end

% Collecting errors with MSE and MAE optimal weights
errors_wmse = y(T2(1)+1:end,:) - sum(a_MSE(1:end-1,:).*Z(T2(1)+1:end,:),2);
errors_wmae = y(T2(1)+1:end,:) - sum(a_MAE(1:end-1,:).*Z(T2(1)+1:end,:),2);
Npred = size(errors_wmse,1);

%% Table 2
% Calculating MSFE and MAFE with MSE and MAE weights.

Weights = {'MSE';'MAE'};
MSE = [sum(errors_wmse.^2,1)/Npred;sum(errors_wmae.^2,1)/Npred];
MAE = [sum(abs(errors_wmse),1)/Npred;sum(abs(errors_wmae),1)/Npred];
Tab = table(Weights, MSE, MAE);

%% Saving weights for Figure 8

tab_empirics = table(Date(T2',1),a_MSE,a_MAE);
writetable(tab_empirics,'Inflation_weights.csv')
%writetable(tab_empirics,'Growth_weights.csv')
%writetable(tab_empirics,'Unemp_weights.csv')
