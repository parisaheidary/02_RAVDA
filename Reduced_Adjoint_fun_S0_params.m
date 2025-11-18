function [J,gradJ]=Reduced_Adjoint_fun_S0_params(S0,S_b,Par0,POD,M_)
% in this function, S0 goes to hydrus to create timeseries, having M_
% and PODs, J and nu_hat and gradJ are calculated. then optimization will
% be performed on this function

pp=1;
vv=2;


z=[0:0.01:0.5,0.52:0.02:1,1.03:0.03:1.51];
lz=length(z);
tt=(1:1:90);


teta_r0=Par0(1,1);
teta_s0=Par0(2,1);
alpha0=Par0(3,1);
x0=Par0(4,1);
m0=Par0(end,1);

if alpha0<0
    alpha0=0.001;
end

if teta_r0<0
    teta_r0=0.0001;
end
if teta_r0>0.1
    teta_r0=0.1;
end

if teta_s0<0.35
    teta_s0=0.35;
end
if teta_s0>0.7
    teta_s0=0.7;
end

% if S0(1,1)>0.25
%     S0(1,1)=0.25;
% end

[t1,S]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0,m0);
% [t1,S]=Hydrus_Run_param_err(S0,teta_r0,teta_s0,alpha0,x0,m0);
load observ.mat observ

[t_d,~,it1]=intersect(tt,t1);
S_obs=observ;
S_est=zeros(1,length(t_d)+1);
% S_est(1,1)=S0(1);
% S_est(1,2:end)=S(1,it1);
S_est(1,1)=mean(S0(1:5,1));
S_est(1,2:end)=mean(S(1:5,it1));



M_1=M_;

lper=length(t_d)+1;

B00=(1/10)*eye(lz,lz); %B0=err^2
R1=1/1250; %R=err^2
aj=1;

J_b=0.5*(S_b-S0)'*B00^-1*(S_b-S0);
% J_b=0;
J=J_b;
for ttt=1:lper
    J=J+0.5*((S_obs(ttt)-S_est(1,ttt))'*R1^-1*(S_obs(ttt)-S_est(1,ttt)));
end


if nargout > 1
    nu_hat_ra=zeros(size(M_1,1),lper);
    
    for t_=lper-1:-1:1
        nu_hat_ra(:,t_) = M_1(:,:,t_+1)' * nu_hat_ra(:,t_+1)+ POD(1,1:vv,pp)'*R1^-1*((S_obs(t_)-S_est(1,t_)));  %becasue H_hat returns P at first row
%         nu_hat_ra(:,t_) = M_1' * nu_hat_ra(:,t_+1)+ POD(1,1:vv,pp)'*R^-1*((S_obs(t_)-S_est(1,t_)));  %becasue H_hat returns P at first row
        
    end
    

%     
%     gradJ=aj*(-B00^-1*(S_b-S0)-POD(:,1:vv,pp)*nu_hat_ra(:,1));

    gradJ=aj*(-POD(:,1:vv,pp)*nu_hat_ra(:,1));



end

fid   =   fopen('costfunS.txt', 'a+');
fprintf(fid, '%6.9f \n', J);
fclose(fid);

% 
figure (2);
% plot(S_b,-z,'b','LineWidth',1.5)
% hold on
plot(S0,-z,'--')
hold on
xlim([0.1 0.45])
% set(gcf,'Position',[950 100 500 550])
set(gcf,'Position',[2500 80 350 350])
% % plot(S_b+0.05,-z,'r')
% hold off

figure(3);
plot(S_obs(1,:),'r','LineWidth',1.5)
hold on
plot(S_est(1,:),'--')
ylim([0 0.45])
set(gcf,'Position',[1500 100 800 300])

% figure(4);
% plot(nu_hat_ra(1,:))
% hold on
% plot(nu_hat_ra(2,:))
% hold off
% set(gcf,'Position',[2300 550 600 200])

% figure(5);
% plot(S_obs-S_est(1,:),'LineWidth',1.5)
% % ylim([-0.02 0.02])
% hold on
% set(gcf,'Position',[1500 450 800 300])

end