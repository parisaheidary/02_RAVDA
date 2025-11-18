function [J,gradJ]=Reduced_Adjoint_fun_S0_params2(S0,Par0,Par_b,POD,M)
% in this function, S0 goes to hydrus to create timeseries, having M_
% and PODs, J and nu_hat and gradJ are calculated. then optimization will
% be performed on this function

pp=1;
vv=6;

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
    teta_r0=0.001;
end
if teta_r0>0.1
    teta_r0=0.1;
end



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



lper=length(t_d)+1;
B0=diag([(1/0.1)*ones(1,1),(1/0.1)*ones(1,1),(1/2)*ones(1,1),(1/0.1)*ones(1,1),(1/0.01)*ones(1,1)]); %B0=err^2
% B0=1/0.9;
R=1/1100;%R=err^2--> err=R^0.5

%%% coef for adjusting/normalizing grad for parameters- parameters are not
%%% in the same units- by experience and trial error
%%%% for wet:
%%%% [1,-0.3,-0.05,-3,-500],
% a_co=[-1,0.27,0.06,3.02,550];
%%%% for dry
%%% [1.8,0.1248,-0.0569,3.02,250]
% a_co=[2.5,0.1248,-0.0569,3.02,250];
a_co=[1.8,0.1248,-0.0569,3.02,250];

a_tr=a_co(1); 
a_ts=a_co(2);
a_a=a_co(3);
a_x=a_co(4);
a_m=a_co(5);

J_b=0.5*(Par_b-Par0)'*B0^-1*(Par_b-Par0);
% J_b=0;
J=J_b;
for ttt=1:lper
    % for ttt=1:lper1
    J=J+0.5*(((S_obs(ttt)-S_est(1,ttt)))'*R^-1*(S_obs(ttt)-S_est(1,ttt)));
    %     J=J+0.5*(((S_obs(t_sel(ttt))-S_est(1,t_sel(ttt))))'*R^-1*(S_obs(t_sel(ttt))-S_est(1,t_sel(ttt))));
end

N=zeros(vv+size(Par0,1),lper);
for ii=1:lper
    N(:,ii)=[POD(1,1:vv,pp)'*R^-1*(S_obs(ii)-S_est(1,ii));zeros(size(Par0,1),1)];
    %     N(:,ii)=[POD(1,1:vv,pp)'*R^-1*(S_obs(t_sel(ii))-S_est(1,t_sel(ii)));zeros(size(Par0,1),1)];
end


if nargout > 1
    nu_hat_ra=zeros(size(M,1),lper);
    
    for t_=lper-1:-1:1
        nu_hat_ra(:,t_) = M(:,:,t_+1)' * nu_hat_ra(:,t_+1)+ N(:,t_+1);  %becasue H_hat returns P at first row
        %         nu_hat_ra(:,t_) = M(:,:,t_+1)' * nu_hat_ra(:,t_+1)+ N(:,t_);
    end
%     
    gradJb=-B0^-1*(Par_b-Par0);
    
%     gradJtr=gradJb(1)-nu_hat_ra(vv+1,1);
    gradJtr=-nu_hat_ra(vv+1,1);
    
%     gradJts=gradJb(2)-nu_hat_ra(vv+2,1);
    gradJts=-nu_hat_ra(vv+2,1);
       
    gradJa=gradJb(3)-nu_hat_ra(vv+3,1);
%         gradJa=-nu_hat_ra(vv+3,1);

    gradJx=gradJb(4)-nu_hat_ra(vv+4,1);
%     gradJx=-nu_hat_ra(vv+4,1);

    gradJm=-nu_hat_ra(end,1);
    

    gradJ=[a_tr*gradJtr;a_ts*gradJts;a_a*gradJa;a_x*gradJx;a_m*gradJm];

       
end
%
fid   =   fopen('costfun.txt', 'a+');
fprintf(fid, '%6.9f \n', J);
fclose(fid);

C=[Par0',gradJ'];
fid   =   fopen('par_grad.txt', 'a+');
fprintf(fid, '%6.9f   %6.9f   %6.9f   %6.9f   \n', C);
fclose(fid);

%
figure(3);
plot(S_obs(1,:),'r','LineWidth',1.5)
hold on
plot(S_est(1,:),'--','LineWidth',1.5)
ylim([0 0.45])
% set(gcf,'Position',[50 100 800 300])
set(gcf,'Position',[2050 100 800 300])
title('Surf')



end