%% GENERATE A RANDOM SIGMA MATRIX
%normaldist = makedist('Normal', 'mu',0,'sigma',1);
%s = random(normaldist,p,p);
%Sigma = s*s';
%eig(Sigma)

%% General paramters
load 'SN_Set2.mat' %or load 'SN_Set2.mat'
p = 5; %number of forecasters.
df = 3;%Degrees of freedom for t-distribution.
R = 10;%5000;%number of replications
N = [20,30];%[20,30,50,100,200,300,400,500,600,700,800,900,1000];%sample size

%Generate Skew Normal Random Variable. Set 1 or Set 2
Lambda = diag(b);
Omega = Sigma + (1-2/pi)*(Lambda*Lambda');

%Set 1 (or 2) of optimal (true) weights 
atheory = Omega\ones(p,1)/(ones(p,1)'/Omega*ones(p,1)); 
display(atheory)

%% Simulations
options = optimset('Display','notify-detailed','LargeScale','off',...
    'MaxFunEvals',5000,'MaxIter',1000);%'TolFun',1e-20,'TolX',1e-20);
asn_MSE = zeros(size(atheory,1),R,numel(N));
asn_MAE = zeros(size(atheory,1),R,numel(N));
at3_MSE = zeros(size(atheory,1),R,numel(N));
at3_MAE = zeros(size(atheory,1),R,numel(N));
amix_MSE = zeros(size(atheory,1),R,numel(N));
amix_MAE = zeros(size(atheory,1),R,numel(N));

for i = 1:numel(N)
    Ni = N(1,i);
    disp(Ni)
    
    parfor r = 1:R
        
        %% Skew Random Forecast Errors
        U = mvnrnd(zeros(p,1),Sigma,Ni);
        tau = mvnrnd(zeros(p,1),eye(p),Ni);
        
        Xi = -sqrt(2/pi)*Lambda*ones(p,1);
        Z = kron(Xi',ones(Ni,1)) + abs(tau)*Lambda+U;
        
        %MSE weights
        Omegahat_SN = (Z'*Z)/Ni;
        asn_MSE(:,r,i) = Omegahat_SN\ones(p,1)/(ones(p,1)'/Omegahat_SN*ones(p,1));
       
        %MAE weights
        asn_init = Omegahat_SN\ones(p,1)/(ones(p,1)'/Omegahat_SN*ones(p,1));
        fsn_mae = @(asn_mae)maeloss(asn_mae,Z);
        asn_mae_raw = fminunc(fsn_mae, asn_init(1:p-1,1),options);
        asn_MAE(:,r,i) = [asn_mae_raw;1-sum(asn_mae_raw)];
        
        %% t_3 forecast errors
        E = mvtrnd(Omega,df,Ni);
        
        %MSE weights
        Omegahat_t3 = (E'*E)/Ni;
        at3_MSE(:,r,i) = Omegahat_t3\ones(p,1)/(ones(p,1)'/Omegahat_t3*ones(p,1));
        
        %MAE weights
        at3_init = Omegahat_t3\ones(p,1)/(ones(p,1)'/Omegahat_t3*ones(p,1));
        ft3_mae = @(at3_mae)maeloss(at3_mae,E);
        at3_mae_raw = fminunc(ft3_mae, at3_init(1:p-1,1),options);
        at3_MAE(:,r,i) = [at3_mae_raw;1-sum(at3_mae_raw)];

    end
    
end

%% Gathering simulated weights and analysis

% Skew normal
amean_SN_diff = squeeze(mean(asn_MSE - asn_MAE,2))';
adev_SN_diff = squeeze(std(asn_MSE-asn_MAE,0,2))';

% T3
amean_t3_diff = squeeze(mean(at3_MSE - at3_MAE,2))';
adev_t3_diff = squeeze(std(at3_MSE-at3_MAE,0,2))';

%% Saving weight differences and standard deviation 
tab_simsnset1 = table(N',amean_SN_diff,adev_SN_diff);
writetable(tab_simsnset1,'sim_snset1.csv')

tab_simt3set1 = table(N',amean_t3_diff,adev_t3_diff);
writetable(tab_simt3set1,'sim_t3set1.csv')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adev_SN_MSE = zeros(numel(N),p-1);
% adev_SN_MAE = zeros(numel(N),p-1);
% amu_SN_MSE = zeros(numel(N),p-1);
% amu_SN_MAE = zeros(numel(N),p-1);
% adev_SN = zeros(numel(N),p-1);
% aplot_SN_diff = zeros(p-1,numel(N),3);
% for k=1:4
%     amu_SN_MSE(:,k) = squeeze(mean(asn_MSE(k,:,:)));
%     amu_SN_MAE(:,k) = squeeze(mean(asn_MAE(k,:,:)));
%     adev_SN_MSE(:,k) = squeeze(std(asn_MSE(k,:,:)));
%     adev_SN_MAE(:,k) = squeeze(std(asn_MAE(k,:,:)));
%     adev_SN(:,k) = squeeze(std(asn_MSE(k,:,:)-asn_MAE(k,:,:)));
%     aplot_SN_diff(k,:,:) = [a_SN_diff(k,:)'-adev_SN(:,k),a_SN_diff(k,:)',a_SN_diff(k,:)'+adev_SN(:,k)];
%     
%     figure;
%     plot(N, squeeze(aplot_SN_diff(k,:,:)))
%     xlim([0 1000]);
%     ylim([-0.1 0.1]);
% end

% adev_t3 = zeros(numel(N),p-1);
% adev_t3_MSE = zeros(numel(N),p-1);
% adev_t3_MAE = zeros(numel(N),p-1);
% amu_t3_MSE = zeros(numel(N),p-1);
% amu_t3_MAE = zeros(numel(N),p-1);
% aplot_t3_diff = zeros(p-1,numel(N),3);
% for k=1:4
%     amu_t3_MSE(:,k) = squeeze(mean(at3_MSE(k,:,:)));
%     amu_t3_MAE(:,k) = squeeze(mean(at3_MAE(k,:,:)));
%     adev_t3(:,k) = squeeze(std(at3_MSE(k,:,:)-at3_MAE(k,:,:)));
%     adev_t3_MSE(:,k) = squeeze(std(at3_MSE(k,:,:)));
%     adev_t3_MAE(:,k) = squeeze(std(at3_MAE(k,:,:)));
%     aplot_t3_diff(k,:,:) = [a_t3_diff(k,:)'-adev_t3(:,k),a_t3_diff(k,:)',a_t3_diff(k,:)'+adev_t3(:,k)];
%     
%     figure;
%     plot(N, squeeze(aplot_t3_diff(k,:,:)))
%     xlim([0 1000]);
%     ylim([-0.1 0.1]);
% end

