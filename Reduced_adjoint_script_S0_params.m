% this code includes the optimization loop- coupled S0 and parameters
clear all;
close all;clc; 
delete('costfun.txt','par_grad.txt')

delete('C:\Hydrus_test\simulation\sim_1\*');
z=[0:0.01:0.5,0.52:0.02:1,1.03:0.03:1.51];
lz=length(z);
tt=(1:1:90);


%% Truth
param_T=[-0.055,-0.165,0.25];
% param_T=[0.055,0.165,0.3738];
St0=(param_T(1)*z'.^2+param_T(2)*(-z')+param_T(3));%+0.03;
teta_rT=0.0883;
teta_sT=0.4662;
alphaT=0.0085;
mT=2.5795;  %Ksat_T=13.19=exp(m);
xT=0.4101;  %nT=1.507=exp(x);

nT=1.507;
Ksat_T=13.19;

[t,St]=Hydrus_Run_param(St0,teta_rT,teta_sT,alphaT,xT,mT);
[t_d,itt,it]=intersect(tt,t);
truth=zeros(lz,length(t_d)+1);
truth(:,1)=St0;
truth(:,2:end)=St(:,it);

% [t_q_T,q_bot_T]=reading_flx_Hy;
% [t_q_T,vtop_T,q_bot_T]=reading_flx_Hy;
% [t_q_T,vtop_T,q_bot_T,inf_T,Ev_T]=reading_flx_Hy;


%% Perturbing obs
a=[-0.06,0.06];  %[-0.04,0.04]
% r = a(1) + 2*a(2)*rand(1,length(tt)); %% in range a
r=(normrnd(0,0.03,[1,90]));  %%%0.04
% r=0;
% observ=truth(1,:);
observ_T=mean(truth(1:5,:));
observ=zeros(1,(length(tt)+1));
observ(1,1)=observ_T(1);
observ(1,2:end)=observ_T(2:end)+r;
idx=find((observ<0.07));
observ(idx)=mean(truth(1:5,idx));
save observ.mat observ

load observ.mat
%%

pp=1;
vv1=2;
vv2=6;
eps1=0.95;
eps2=0.95;
%%% perturbed using the range of each param
del=1/100;
dtr=del*(0.1-0);
dts=del*(0.7-0.2);
da=del*(0.06-0);
dx=del*(0.8329-0);  
dm=del*(144-12);

num_par=5;
teta_r0=0.008; 
teta_s0=0.389;  
teta_s0_1=0.6; 
alpha0=0.023;
x0=log(1.23);% n=exp(x)
m0=log(2.7183);  %log(13.19)=2.5795    K=exp(m)   log(13.19)=2.5795


S_0=St0(1,1)*ones(lz,1);
% S_0=0.38*ones(lz,1);

%% OL
% to have S_mean
[t0,S00]=Hydrus_Run_param(S_0,teta_r0,teta_s0,alpha0,x0,m0);
% [t0,S00]=Hydrus_Run_param_err(S_0,teta_r0,teta_s0,alpha0,x0,m0);
[t_d0,~,it0]=intersect(tt,t0);
theta0=zeros(lz,length(t_d0)+1);
% theta0(:,1)=S_0;
% theta0(:,2:end)=S00(:,it0);
theta0(:,1)=S_0;
theta0(:,2:end)=S00(:,it0);
surf_OL=mean(theta0(1:5,:));

% [t_q_ol,vtop_ol,q_bot_ol,inf_ol,Ev_ol]=reading_flx_Hy;

% [t_q_ol,vtop_ol,q_bot_ol,inf_ol,Ev_ol]=reading_flx_Hy_err;


S0=S_0;

%% loaded P and Ms that I created prior to save time-
% because I am not updating them at each iteration

filename1=strcat('P_M_ttaxm_S_',num2str(eps1),'_','r',num2str(vv1),'p',num2str(pp),'_',num2str(S0(1)),'_tr0_',num2str(teta_r0),'_ts0_',num2str(teta_s0_1),'_a0_',num2str(alpha0),'_x0_',num2str(x0),'_m0_',num2str(m0_1),'.mat');
fn1=fullfile('C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\PM_matfiles',filename1);
load(fn1,'POD','M_')

filename2=strcat('P_M_ttaxm_',num2str(eps2),'_dtr_',num2str(dtr),'_dts_',num2str(dts),'_da_',num2str(da),'_dx_',num2str(dx),'_','_dm_',num2str(dm),'_','r',num2str(vv2),'p',num2str(pp),'_',num2str(S0(1)),'_tr0_',num2str(teta_r0),'_ts0_',num2str(teta_s0),'_a0_',num2str(alpha0),'_x0_',num2str(x0),'_m0_',num2str(m0),'.mat');
fn2=fullfile('C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\PM_matfiles',filename2);
load(fn2,'POD_comb','M')





options=optimset('Display','iter','MaxIter',1000,'TolFun',0.005, 'MaxFunEvals', 100,'TolX', 0.005, 'LargeScale','off','GradObj','on','DerivativeCheck','off','HessUpdate','bfgs','FinDiffType','forward');%,'LineSearchType','quadcubic');%,'TolGradCon', 1e-6, 'TolCon', 1e-6);%,'Diagnostics','on');%,'PlotFcns',@myplotfun1);'CheckGradients','on'  %%'HessUpdate','dfp' 'TolGradCon', 1e-10

S_b(:,1)=S0;

% Par0=[0.05;0.45;0.03;0.42;3.73]; % ave of ranges
Par0=[teta_r0;teta_s0;alpha0;x0;m0]; 
% Par0=[teta_rT;teta_sT;alphaT;xT;mT]; 
Par_b(:,1)=Par0;
Par_b0(:,1)=Par0;

n=0;
J_i1(1)=100;J_i2(1)=100;miu=100;

while miu>1e-25
        
    n=n+1
    
    [S0_up,fval,~,~,grad,hessian]=fminunc(@(S0)Reduced_Adjoint_fun_S0_params(S0,S_b(:,n),Par0,POD,M_),S0,options);
    
    S_b0(:,n+1)=S0_up;
    S0=S0_up;
    S_b(:,n+1)=S0;
    gradS(:,n)=grad;
    hessian_s(:,:,n)=hessian;
    J_i1(n+1)=fval;
    miu1(n+1)=abs((J_i1(n+1)-J_i1(n)))/max(abs(J_i1(n+1)),1);
    
    [Par0_up,fval,exitflag,output,grad,hessian]=fminunc(@(Par0)Reduced_Adjoint_fun_S0_params2(S0,Par0,Par_b(:,n),POD_comb,M),Par0,options);
    
    J_i2(n+1)=fval;
    miu2(n+1)=abs((J_i2(n+1)-J_i2(n)))/max(abs(J_i2(n+1)),1);
    
    %%
    if Par0_up(3,1)<0
        Par0_up(3,1)=0.001;
    end
    
    if Par0_up(1,1)<0
        Par0_up(1,1)=0.0001;
    end
    if Par0_up(1,1)>0.1
        Par0_up(1,1)=0.1;
    end
    
    if Par0_up(2,1)<0.2
        Par0_up(2,1)=0.2;
    end
    if Par0_up(2,1)>0.7
        Par0_up(2,1)=0.7;
    end
%%        
    Par_b0(:,n+1)=Par0_up;
    Par0=Par0_up;
    Par_b(:,n+1)=Par0;
    grad_par(:,n)=grad;
    hessian_par(:,:,n)=hessian;
    
    miu(n+1)=max(abs(miu1(n+1)),abs(miu2(n+1)));
    miu(n+1)

end
Par0=Par0_up;
S0_est=S0;
teta_r_est=Par0(end-4,1);
teta_s_est=Par0(end-3,1);
alpha_est=Par0(end-2,1);
n_est=exp(Par0(end-1,1));
K_est=exp(Par0(end,1));



[t_est,S_est0]=Hydrus_Run_param(S0_est,teta_r_est,teta_s_est,alpha_est,Par0(end-1,1),Par0(end,1));
% [t_est,S_est0]=Hydrus_Run_param_err(S0_est,teta_r_est,teta_s_est,alpha_est,Par0(end-1,1),Par0(end,1));
[t_d00,~,it_est]=intersect(tt,t_est);
S_est=zeros(lz,length(t_d0)+1);
S_est(:,1)=S0_est;
S_est(:,2:end)=S_est0(:,it_est);
S_est_surf=mean(S_est(1:5,:));

% [t_q_est,q_bot_est]=reading_flx_Hy;
% [t_q_est,q_bot_est]=reading_flx_Hy_err;


figure
plot(St0,-z,'r','LineWidth', 2)%8
hold on
% plot(0.35*ones(93,1),-z,'b','LineWidth', 2)
plot(S_b(:,1),-z,'b','LineWidth', 2)
hold on
plot(S_b(:,n+1),-z,'k--','LineWidth', 2)
% plot(S0_up,-z,'k--','LineWidth',2)
% hold off
% legend({'True Initial soil moisture profile',['RA-estimated' newline 'initial soil moisture profile'],'initial guess','S-b'},'fontsize',10)
legend({'True SM','Guessed SM','Estimated SM'},'fontweight','bold','fontsize',12)%, 53, 'Orientation','horizontal')
xlabel('Initial SM Profile','fontweight','bold','fontsize',12);%55
% xlim([0.25 0.42])
set(gca,'TickDir','out')
ylabel(' Depth (m)','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4
xlim([0.1 0.45]);
% set(gcf,'Position',[-30 -30 1150 2050])
% set(gcf,'Position',[50 50 650 800])
set(gcf,'Position',[100 100 550 600])
legend boxoff


figure
plot(Par_b(1,:),'--k','LineWidth', 2);
hold on
plot(teta_rT*ones(1,n+1),'LineWidth', 2);
hold off
% legend({'Estimated Ksat','True Ksat'},'fontweight','bold','fontsize',12);%, 53, 'Orientation','horizontal')
legend({'Estimated teta_r','True teta_r'},'fontweight','bold','fontsize',12)
xlabel('Iterations','fontweight','bold','fontsize',12);%55
% ylim([12 21])
xlim([1,n+1])
set(gca,'TickDir','out')
% ylabel(' Ksat cm/day','fontweight','bold','fontsize',12);
ylabel(' teta_r cm/cm','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4


figure
plot(Par_b(2,:),'--k','LineWidth', 2);
hold on
plot(teta_sT*ones(1,n+1),'LineWidth', 2);
hold off
% legend({'Estimated Ksat','True Ksat'},'fontweight','bold','fontsize',12);%, 53, 'Orientation','horizontal')
legend({'Estimated teta_s','True teta_s'},'fontweight','bold','fontsize',12)
xlabel('Iterations','fontweight','bold','fontsize',12);%55
% ylim([12 21])
xlim([1,n+1])
set(gca,'TickDir','out')
% ylabel(' Ksat cm/day','fontweight','bold','fontsize',12);
ylabel(' teta_s cm/cm','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4


figure
plot(Par_b(3,:),'--k','LineWidth', 2);
hold on
% plot(Ksat_T*ones(1,n+1),'LineWidth', 2);
plot(alphaT*ones(1,n+1),'LineWidth', 2);
hold off
% legend({'Estimated Ksat','True Ksat'},'fontweight','bold','fontsize',12);%, 53, 'Orientation','horizontal')
legend({'Estimated alpha','True alpha'},'fontweight','bold','fontsize',12)
xlabel('Iterations','fontweight','bold','fontsize',12);%55
% ylim([12 21])
xlim([1,n+1])
set(gca,'TickDir','out')
% ylabel(' Ksat cm/day','fontweight','bold','fontsize',12);
ylabel(' alpha 1/cm','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4

figure
plot(exp(Par_b(end,:)),'--k','LineWidth', 2);
hold on
plot(Ksat_T*ones(1,n+1),'LineWidth', 2);
% plot(alphaT*ones(1,n+1),'LineWidth', 2);
hold off
legend({'Estimated Ksat','True Ksat'},'fontweight','bold','fontsize',12);%, 53, 'Orientation','horizontal')
% legend({'Estimated alpha','True alpha'},'fontweight','bold','fontsize',12)
xlabel('Iterations','fontweight','bold','fontsize',12);%55
% ylim([12 20])
xlim([1,n+1])
set(gca,'TickDir','out')
ylabel(' Ksat cm/day','fontweight','bold','fontsize',12);
% ylabel(' alpha 1/cm','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4

figure
plot(exp(Par_b(end-1,:)),'--k','LineWidth', 2);
hold on
plot(nT*ones(1,n+1),'LineWidth', 2);
% plot(alphaT*ones(1,n+1),'LineWidth', 2);
hold off
legend({'Estimated n','True n'},'fontweight','bold','fontsize',12);%, 53, 'Orientation','horizontal')
% legend({'Estimated alpha','True alpha'},'fontweight','bold','fontsize',12)
xlabel('Iterations','fontweight','bold','fontsize',12);%55
% ylim([12 17])
xlim([1,n+1])
set(gca,'TickDir','out')
ylabel(' n [-]','fontweight','bold','fontsize',12);
% ylabel(' alpha 1/cm','fontweight','bold','fontsize',12);
set(gca,'FontName','Arial','FontSize',12,'FontWeight','Bold','LineWidth', 2);%4




figure
% plot(truth(1,:),'-.g*')%,'LineWidth',1.5)
plot(observ_T(1,:),'-.g*')%,'LineWidth',1.5)
% plot(y_t(1,:),'-.g*')%,'LineWidth',1.5)
% hold on
% plot(observ(1,:),'-.r*')%,'LineWidth',1.5)
hold on
% plot(theta0(1,:),'b','LineWidth',1.5)
plot(surf_OL(1,:),'b','LineWidth',1.5)
hold on
plot(S_est_surf(1,:),'k--','LineWidth',2)
% legend({'True','obs','OL','RA-VDA'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
legend({'True','OL','RA-VDA'},'fontweight','bold','fontsize',20,'Orientation','horizontal')
ylabel('SM at the surface')
xlabel('Day')
ylim([0.05 0.4])
set(gcf,'position',[100,200,1110,500])
set(gca,'FontName','Times New Roman','fontweight','bold','FontSize',20,'LineWidth', 2);
legend boxoff


figure
plot(truth(51,:),'-.g*')%,'LineWidth',1.5)
% plot(truth(51,:),'-.r*')%,'LineWidth',1.5)
hold on
% plot(theta0(1,:),'b','LineWidth',1.5)
plot(theta0(51,:),'b','LineWidth',1.5)
hold on
plot(S_est(51,:),'k--','LineWidth',2)
legend({'True','OL','RA-VDA'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
ylabel('SM at the Middle')
xlabel('Day')
ylim([0.05 0.4])
set(gcf,'position',[100,200,1110,500])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
legend boxoff


figure
plot(truth(93,:),'-.g*')%,'LineWidth',1.5)
% plot(truth(93,:),'-.r*')%,'LineWidth',1.5)
hold on
% plot(theta0(1,:),'b','LineWidth',1.5)
plot(theta0(93,:),'b','LineWidth',1.5)
hold on
plot(S_est(93,:),'k--','LineWidth',2)
legend({'True','OL','RA-VDA'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
ylabel('SM at the Bottom')
xlabel('Day')
ylim([0.05 0.4])
set(gcf,'position',[100,200,1110,500])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
legend boxoff

figure;
plot(t_q_T,q_bot_T,'g','LineWidth',2)
hold on
plot(t_q_ol,q_bot_ol,'b','LineWidth',2)
hold on
plot(t_q_est,q_bot_est,'k--','LineWidth',2)
legend({'True','OL','RA-VDA'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
ylabel('Bottom flux')
xlabel('Day')
% ylim([0.05 0.4])
set(gcf,'position',[100,200,1110,500])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
legend boxoff

[t_q_io,i_qTo,i_qo]=intersect(t_q_T,t_q_ol);
[t_q_i,i_qT,i_qe]=intersect(t_q_T,t_q_est);
figure;
scatter(q_bot_T(i_qT),q_bot_est(i_qe),'k*');
ylabel('Estimated flux')
xlabel('True flux')
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
% xlim([-0.2 0])
% ylim([-0.2 0])

RMSE_qo=sqrt(mean(q_bot_T(i_qTo)-q_bot_ol(i_qo)).^2)
RMSE_qe=sqrt(mean(q_bot_T(i_qT)-q_bot_est(i_qe)).^2)
cor_q=corrcoef(q_bot_T(i_qT),q_bot_est(i_qe))


tr_ol=(Par_b(1,1)-teta_rT)*100/teta_rT;
tr_est=(teta_r_est-teta_rT)*100/teta_rT;

ts_ol=(Par_b(2,1)-teta_sT)*100/teta_sT;
ts_est=(teta_s_est-teta_sT)*100/teta_sT;

a_ol=(Par_b(3,1)-alphaT)*100/alphaT;
a_est=(alpha_est-alphaT)*100/alphaT;

% x_ol=(Par_b(4,1)-xT)*100/xT;
% x_est=(Par0(end-1,1)-xT)*100/xT;
% 
% m_ol=(Par_b(5,1)-mT)*100/mT;
% m_est=(Par0(end,1)-mT)*100/mT;

n_ol=(exp(Par_b(4,1))-nT)*100/nT;
n_es=(n_est-nT)*100/nT;

k_ol=(exp(Par_b(5,1))-Ksat_T)*100/Ksat_T;
k_est=(K_est-Ksat_T)*100/Ksat_T;

RMSE_S0_ol=sqrt(mean((S_b(:,1)-St0(:,1)).^2))
RMSE_S0_es=sqrt(mean((S_b(:,end)-St0(:,1)).^2))

RMSE_ol=sqrt(mean((surf_OL(1,:)-observ_T(1,:)).^2))
RMSE_es=sqrt(mean((S_est_surf(1,:)-observ_T(1,:)).^2))

% RMSE_ol=sqrt(mean((surf_OL(1,:)-observ(1,:)).^2))
% RMSE_es=sqrt(mean((S_est_surf(1,:)-observ(1,:)).^2))

RMSE_ol1=sqrt(mean((surf_OL(1,1:60)-observ_T(1,1:60)).^2))
RMSE_ol2=sqrt(mean((surf_OL(1,61:90)-observ_T(1,61:90)).^2))

RMSE_es1=sqrt(mean((S_est_surf(1,1:60)-observ_T(1,1:60)).^2))
RMSE_es2=sqrt(mean((S_est_surf(1,61:90)-observ_T(1,61:90)).^2))

JJ1=J_i1/J_i1(1);
JJ2=J_i2/J_i2(1);
plot(JJ1,'k-*','LineWidth',1.5)
hold on
plot(JJ2,'b-+','LineWidth',1.5)
legend('S0','parameters')
ylabel('normalized J')
xlabel('outer iterations')
ylabel('normalized J, J/J0')

