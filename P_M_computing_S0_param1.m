function [POD,M_]=P_M_computing_S0_param1(S0,teta_r0,teta_s0,alpha0,x0,m0)

tic
z=[0:0.01:0.5,0.52:0.02:1,1.03:0.03:1.51];
lz=length(z);
num_node=lz;
tt=(1:1:90);


%%%name='P_M_computing_S0_param1';

pp=1;
r=(1:2);
v=length(r);
eps=0.95;


delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_M,theta]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0,m0);
[t_d,~,it]=intersect(tt,t_M);
theta_d=zeros(lz,length(t_d)+1);
theta_d(:,1)=S0;
theta_d(:,2:end)=theta(:,it);


S_bar=zeros(lz,1);
POD=zeros(lz,6,v);

% for iii=1:1  %% end is in line 96
iii=1; %%% every day snapshot


S_bar(:,iii)=mean(theta_d(:,1:iii:end),2);
E_S=theta_d(:,1:iii:end)-S_bar(:,iii);


co_var_S=E_S'*E_S;
[Z_S,D_S]=eigs(co_var_S);

%% for wet
% %
% if Z_S(1,1)*Z_S(1,2)>0 && Z_S(1,1)>0
%    Z_S(:,1)=-Z_S(:,1);
%    Z_S(:,2)=-Z_S(:,2);
% end
%
% if Z_S(1,1)*Z_S(1,2)<0
%     if Z_S(1,1)>0
%     Z_S(:,1)=-Z_S(:,1);
%     else
%    Z_S(:,2)=-Z_S(:,2);
%     end
% end

%% for dry
if Z_S(1,1)*Z_S(1,2)>0 && Z_S(1,1)<0
    Z_S(:,1)=-Z_S(:,1);
    Z_S(:,2)=-Z_S(:,2);
end

if Z_S(1,1)*Z_S(1,2)<0
    if  Z_S(1,2)>0
        Z_S(:,1)=-Z_S(:,1);
    else
        Z_S(:,2)=-Z_S(:,2);
    end
end

%%


POD(:,:,iii)=(E_S(:,:)*Z_S(:,:))*(abs(D_S))^-0.5;

%     filename11=strcat(path,'\',num2str(iii));
%     save(filename11,'Z_S','D_S','co_var_S','E_S');

% end


%% M_bar w/ respect to S

S_pert=zeros(num_node,length(t_d),v);
for m=1:v
    for mm=1:length(t_d)
        S_pert(:,mm,m)=(theta_d(:,mm)+eps*POD(:,m,pp));
        S_pert(:,mm,m)=abs(S_pert(:,mm,m));

    end
end
%%%% for stability
S_pert(S_pert<0.089)=0.089+0.03;
S_pert(S_pert>0.4702)=0.4702;


S_prime=zeros(lz,length(t_d)+1,v);
for m=1:v
    for mm=1:length(t_d)
        delete('C:\Hydrus_test\simulation\sim_param_1\*');
        [t1,SS]=Hydrus_Run_param(S_pert(:,mm,m),teta_r0,teta_s0,alpha0,x0,m0); %t1 and whole SS is not important, becasue we need just SS at time 2
        [~,~,it1]=intersect(tt,t1);
        S_prime(:,mm+1,m)=SS(:,it1(1)); %it1(1) because Hydrus starts from next time step and not initial one but, sometimes need to be it1(2)


        %         S_prime(:,mm+1,m)=SS(:,2);

    end
end

M_=zeros(v,v,length(t_d));ratio=zeros(lz,v,length(t_d)+1);
for mmm=2:(length(t_d)+1)
    for vv=1:v

        ratio(:,vv,mmm)=((S_prime(:,mmm,vv)-theta_d(:,mmm))/eps);

    end
    M_(:,:,mmm)=POD(:,1:vv,1)'*ratio(:,1:vv,mmm);
end


toc

%% plots

con_num=zeros(1,length(t_d)+1);
% e=zeros(2,91);
% con_num1=zeros(1,length(t_dk)+1);
for x=2:91
    e(:,x)=real(eig(M_(:,:,x)));
    %     variance(:,i)=1./eigen(:,i);
    con_num(x)=cond(M_(:,:,x));
    con_num1(x)=cond(M_(:,:,x),1);
    rcon_num(x)=rcond(M_(:,:,x));
    iM(:,:,x)=inv(M_(:,:,x));
    %     con_num1(i)=cond(M(:,:,i));
end

figure;plot(con_num(2:91),'LineWidth',2)
figure;plot(e(1,2:91));hold on;plot(e(2,2:91));hold off;



figure;
plot(POD(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD(:,3,iii),-z,'--','LineWidth',2)
hold on
% plot(POD(:,1,iii)+POD(:,2,iii)+POD(:,3,iii),-z,'LineWidth',2)
plot(POD(:,1,iii)+POD(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD1','POD2','comb')

rel_eng1=zeros(size(D_S,1),1);
for m=1:size(D_S,1)
    rel_eng1(m)=D_S(m,m)/trace(D_S);
end

sai1(1)=rel_eng1(1);
for m=2:size(D_S,1)
    sai1(m)=rel_eng1(m)+sai1(m-1);

end


figure
plot(100*sai1,'-ok','LineWidth', 3);
legend({'POD'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
xlabel('No. of POD modes','fontweight','bold','fontsize',18)
ylabel('%','fontweight','bold','fontsize',18)
ylim([90 102])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
% legend boxoff

filename=strcat('P_M_ttaxm_S_',num2str(eps),'_','r',num2str(vv),'p',num2str(pp),'_',num2str(S0(1)),'_tr0_',num2str(teta_r0),'_ts0_',num2str(teta_s0),'_a0_',num2str(alpha0),'_x0_',num2str(x0),'_m0_',num2str(m0),'.mat');
fn=fullfile('C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\PM_matfiles',filename);
save(fn);

end


