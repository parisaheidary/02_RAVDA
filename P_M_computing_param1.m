% function [POD_k,M_,M_k,M]=P_M_computing_param(S0,alpha0,K0)
% function [POD_k,M_,M_k,M]=P_M_computing_param(S0,alpha0,n0,K0)
% function [POD_comb,M_,M_a,M_m,M]=P_M_computing_param(S0,alpha0,n0,m0)
function [POD_comb,M_,M_a,M_m,M]=P_M_computing_param1(S0,teta_r0,teta_s0,alpha0,x0,m0)

name='P_M_computing_param1';

tic
z=[0:0.01:0.5,0.52:0.02:1,1.03:0.03:1.51];
lz=length(z);
num_node=lz;
tt=(1:1:90);

num_par=5;

del=1/100;
pp=1;
r=(1:6);
v=length(r);
eps=0.8;
dtr=del*(0.1-0);
dts=del*(0.7-0.2);
da=del*(0.06-0);
% dn=0.001;
dx=del*(0.8329-0);  
dm=del*(log(144)-log(12));


delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_M,theta]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0,m0);
[t_d,~,it]=intersect(tt,t_M);
theta_d=zeros(lz,length(t_d)+1);
theta_d(:,1)=S0;
theta_d(:,2:end)=theta(:,it);

% work_dir='C:\Users\parisaheidary\Desktop\testData';
% path=strcat(work_dir);

%% POD w/ teta_r

tetar0_pert=teta_r0+dtr;
% al_pert=al+da;
delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_Mtr,theta_tr]=Hydrus_Run_param(S0,tetar0_pert,teta_s0,alpha0,x0,m0);
[t_dtr,~,ittr]=intersect(tt,t_Mtr);
theta_dtr=zeros(lz,length(t_dtr)+1);
theta_dtr(:,1)=S0;
theta_dtr(:,2:end)=theta_tr(:,ittr);


iii=1;
if teta_r0>tetar0_pert
    E_tr=theta_d(:,1:iii:end)-theta_dtr(:,1:iii:end);
else
    E_tr=theta_dtr(:,1:iii:end)-theta_d(:,1:iii:end);
end


co_var_tr=E_tr'*E_tr;
[Z_tr,D_tr]=eigs(co_var_tr,10);

POD_tr(:,:,iii)=(E_tr(:,:)*Z_tr(:,:))*(abs(D_tr))^-0.5;


%% POD w/ teta_s

tetas0_pert=teta_s0+dts;
% al_pert=al+da;
delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_Mts,theta_ts]=Hydrus_Run_param(S0,teta_r0,tetas0_pert,alpha0,x0,m0);
[t_dts,~,itts]=intersect(tt,t_Mts);
theta_dts=zeros(lz,length(t_dts)+1);
theta_dts(:,1)=S0;
theta_dts(:,2:end)=theta_ts(:,itts);


iii=1;
if teta_s0>tetas0_pert
    E_ts=theta_d(:,1:iii:end)-theta_dts(:,1:iii:end);
else
    E_ts=theta_dts(:,1:iii:end)-theta_d(:,1:iii:end);
end


co_var_ts=E_ts'*E_ts;
[Z_ts,D_ts]=eigs(co_var_ts,10);

POD_ts(:,:,iii)=(E_ts(:,:)*Z_ts(:,:))*(abs(D_ts))^-0.5;
%

%% POD w/ alpha

alpha0_pert=alpha0+da;
% al_pert=al+da;
delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_Ma,theta_a]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0_pert,x0,m0);
[t_da,~,ita]=intersect(tt,t_Ma);
theta_da=zeros(lz,length(t_da)+1);
theta_da(:,1)=S0;
theta_da(:,2:end)=theta_a(:,ita);

iii=1;
if alpha0>alpha0_pert
    E_a=theta_d(:,1:iii:end)-theta_da(:,1:iii:end);
else
    E_a=theta_da(:,1:iii:end)-theta_d(:,1:iii:end);
end

co_var_a=E_a'*E_a;
[Z_a,D_a]=eigs(co_var_a,10);


POD_a(:,:,iii)=(E_a(:,:)*Z_a(:,:))*(abs(D_a))^-0.5;

% end


%% POD w/ x

x0_pert=x0+dx;
delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_Mx,theta_x]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0_pert,m0);
[t_dx,~,itx]=intersect(tt,t_Mx);
theta_dx=zeros(lz,length(t_dx)+1);
theta_dx(:,1)=S0;
theta_dx(:,2:end)=theta_x(:,itx);
% 
iii=1;
if x0>x0_pert
    E_x=theta_d(:,1:iii:end)-theta_dx(:,1:iii:end);
else
    E_x=theta_dx(:,1:iii:end)-theta_d(:,1:iii:end);
end

co_var_x=E_x'*E_x;
[Z_x,D_x]=eigs(co_var_x,10);
POD_x(:,:,iii)=(E_x(:,:)*Z_x(:,:))*(abs(D_x))^-0.5;
% 
% % % end

%% POD w/ m

m0_pert=m0+dm;
delete('C:\Hydrus_test\simulation\sim_param_1\*');
[t_Mm,theta_m]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0,m0_pert);
[t_dm,~,itm]=intersect(tt,t_Mm);
theta_dm=zeros(lz,length(t_dm)+1);
theta_dm(:,1)=S0;
theta_dm(:,2:end)=theta_m(:,itm);


iii=1;
if m0_pert<m0  
    E_m=theta_d(:,1:iii:end)-theta_dm(:,1:iii:end);
else
    E_m=theta_dm(:,1:iii:end)-theta_d(:,1:iii:end);
end

co_var_m=E_m'*E_m;
[Z_m,D_m]=eigs(co_var_m,10);

POD_m(:,:,iii)=(E_m(:,:)*Z_m(:,:))*(abs(D_m))^-0.5;

% end

%% POD w/ combination of PODs

E_tr1=E_tr-mean(E_tr);
E_tr2=E_tr-mean(E_tr,2);

E_ts1=E_ts-mean(E_ts);
E_ts2=E_ts-mean(E_ts,2);

E_a1=E_a-mean(E_a);
E_a2=E_a-mean(E_a,2);

E_x1=E_x-mean(E_x);
E_x2=E_x-mean(E_x,2);

E_m1=E_m-mean(E_m);
E_m2=E_m-mean(E_m,2);

E_comb=[E_tr1(:,2:end),E_ts1(:,2:end),E_a1(:,2:end),E_x1(:,2:end),E_m1(:,2:end)];


co_var_comb=E_comb'*E_comb;
[Z_comb,D_comb]=eigs(co_var_comb,10);


POD_comb(:,:,iii)=(E_comb(:,:)*Z_comb(:,:))*(abs(D_comb))^-0.5;


%% M_tetar
M_tr=zeros(v,length(t_d));
% ratio_tr=E_tr1/dtr;
ratio_tr=E_tr2/dtr;
% M_tr=POD_tr(:,1:v,1)'*ratio_tr;
M_tr=POD_comb(:,1:v,1)'*ratio_tr;

%% M_tetas
M_ts=zeros(v,length(t_d));
% ratio_ts=E_ts1/dts;
ratio_ts=E_ts2/dts;
% M_ts=POD_ts(:,1:v,1)'*ratio_ts;
M_ts=POD_comb(:,1:v,1)'*ratio_ts;

%% M_a
ratio_a=E_a2/da;
% M_a=POD_a(:,1:v,1)'*ratio_a;
M_a=POD_comb(:,1:v,1)'*ratio_a;

%% M_x

% ratio_x=E_x1/dx;
ratio_x=E_x2/dx;
% M_x=POD_x(:,1:v,1)'*ratio_x;
M_x=POD_comb(:,1:v,1)'*ratio_x;

%% M_m
ratio_m=E_m2/dm;
% M_m=POD_m(:,1:v,1)'*ratio_m;
M_m=POD_comb(:,1:v,1)'*ratio_m;


%% M_bar w/ respect to S

S_pert=zeros(num_node,length(t_d),v);
for m=1:v
    for mm=1:length(t_d)
        S_pert(:,mm,m)=(theta_d(:,mm)+eps*POD_comb(:,m,pp));
        S_pert(:,mm,m)=abs(S_pert(:,mm,m));
        
    end
end

%%% for stabilty of the solution
S_pert(S_pert<(teta_r0+0.005))=(teta_r0+0.005)+0.03;
S_pert(S_pert>(teta_s0+0.005))=teta_s0+0.005;


S_prime=zeros(lz,length(t_d)+1,v);
for m=1:v
    for mm=1:length(t_d)
        delete('C:\Hydrus_test\simulation\sim_sk_1\*')
%t1 and whole SS is not important, becasue we need just SS at time 2
%         [t1,SS]=Hydrus_Run_param(S_pert(:,mm,m),alpha0,n0,m0);
        [t1,SS]=Hydrus_Run_param(S_pert(:,mm,m),teta_r0,teta_s0,alpha0,x0,m0);
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
     M_(:,:,mmm)=POD_comb(:,1:vv,1)'*ratio(:,1:vv,mmm);

end

%% M

M=zeros(v+num_par,v+num_par,length(t_d)+1);
% M=zeros(v+num_par,v+num_par,length(t_sel));
for i=2:(length(t_d)+1)  
        M(:,:,i)=[M_(:,:,i),M_tr(:,i),M_ts(:,i),M_a(:,i),M_x(:,i),M_m(:,i);zeros(num_par,v),eye(num_par)];
end

%% plots

con_num=zeros(1,length(t_d)+1);
e=zeros(v,91);
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
% figure;plot(con_num1(2:91),'LineWidth',2)

toc
% for iii=1:1
iii=1;

figure;
plot(POD_tr(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD_tr(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_tr(:,3,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_tr(:,1,iii)+POD_tr(:,2,iii)+POD_tr(:,3,iii),-z,'LineWidth',2)
plot(POD_tr(:,1,iii)+POD_tr(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD_{tr1}','POD_{tr2}','comb')

figure;
plot(POD_ts(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD_ts(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_ts(:,3,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_ts(:,1,iii)+POD_ts(:,2,iii)+POD_ts(:,3,iii),-z,'LineWidth',2)
plot(POD_ts(:,1,iii)+POD_ts(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD_{ts1}','POD_{ts2}','comb')

figure;
plot(POD_a(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD_a(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_a(:,3,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_a(:,1,iii)+POD_a(:,2,iii)+POD_a(:,3,iii),-z,'LineWidth',2)
plot(POD_a(:,1,iii)+POD_a(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD_a1','POD_a2','comb')

figure;
plot(POD_x(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD_x(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_n(:,3,iii),-z,'--','LineWidth',2)
% hold on
% plot(POD_n(:,1,iii)+POD_n(:,2,iii)+POD_n(:,3,iii),-z,'LineWidth',2)
plot(POD_x(:,1,iii)+POD_x(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD_x1','POD_x2','comb_x')

figure;
plot(POD_m(:,1,iii),-z,'--','LineWidth',2)
hold on
plot(POD_m(:,2,iii),-z,'--','LineWidth',2)
hold on
% plot(POD_m(:,3,iii),-z,'--','LineWidth',2)
% hold on
% plot(POD_m(:,1,iii)+POD_m(:,2,iii)+POD_m(:,3,iii),-z,'LineWidth',2)
plot(POD_m(:,1,iii)+POD_m(:,2,iii),-z,'LineWidth',2)
hold off
% legend('POD1','POD2','POD3','comb')
legend('POD_m1','POD_m2','comb_m')

figure;
plot(POD_comb(:,1),-z,'--','LineWidth',2)
hold on
plot(POD_comb(:,2),-z,'--','LineWidth',2)
hold on
% plot(POD_comb(:,3),-z,'--','LineWidth',2)
hold on
% plot(POD_comb(:,1)+POD_comb(:,2)+POD_comb(:,3),-z,'LineWidth',2)
plot(POD_comb(:,1)+POD_comb(:,2),-z,'LineWidth',2)
hold off
% legend('POD_comb1','POD_comb2','POD_comb3','comb')
legend('POD_comb1','POD_comb2','comb')
%%% plot(POD_comb(:,1)+POD_comb(:,2)+POD_comb(:,3)+POD_comb(:,4)+POD_comb(:,5)+POD_comb(:,6)+POD_comb(:,7),-z,'LineWidth',2)


rel_eng1=zeros(size(D_a,1),1);
rel_eng2=zeros(size(D_x,1),1);
rel_eng4=zeros(size(D_comb,1),1);
rel_eng3=zeros(size(D_m,1),1);
rel_eng5=zeros(size(D_tr,1),1);
rel_eng6=zeros(size(D_ts,1),1);
 for m=1:size(D_a,1)
  rel_eng1(m)=D_a(m,m)/trace(D_a);
  rel_eng2(m)=D_x(m,m)/trace(D_x);
  rel_eng3(m)=D_m(m,m)/trace(D_m);
  rel_eng4(m)=D_comb(m,m)/trace(D_comb);
  rel_eng5(m)=D_tr(m,m)/trace(D_tr);
  rel_eng6(m)=D_ts(m,m)/trace(D_ts);
 end
 
 
sai1(1)=rel_eng1(1);
sai2(1)=rel_eng2(1);
sai3(1)=rel_eng3(1);
sai4(1)=rel_eng4(1);
sai5(1)=rel_eng5(1);sai6(1)=rel_eng6(1);
for m=2:size(D_a,1)
    sai1(m)=rel_eng1(m)+sai1(m-1);
    sai2(m)=rel_eng2(m)+sai2(m-1);
    sai3(m)=rel_eng3(m)+sai3(m-1);
    sai4(m)=rel_eng4(m)+sai4(m-1);
    sai5(m)=rel_eng5(m)+sai5(m-1);
    sai6(m)=rel_eng6(m)+sai6(m-1);
end
figure
plot(100*sai5,'-ok','LineWidth', 3);
hold on
plot(100*sai6,'-oc','LineWidth', 3);
hold on
plot(100*sai1,'-or','LineWidth', 3);
hold on
plot(100*sai2,'-sM','LineWidth', 3);
hold on
plot(100*sai3,'-sb','LineWidth', 3);
hold on
plot(100*sai4,'-*g','LineWidth', 3);
legend({'POD_{tr}','POD_{ts}','POD_{a}','POD_{x}','POD_{m}','POD_{comb}'},'fontweight','bold','fontsize',18,'Orientation','horizontal')
xlabel('No. of POD modes','fontweight','bold','fontsize',18)
ylabel('%','fontweight','bold','fontsize',18)
ylim([90 102])
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth', 1.2);
% legend boxoff


%% Save Files
filename=strcat('P_M_ttaxm_',num2str(eps),'_dtr_',num2str(dtr),'_dts_',num2str(dts),'_da_',num2str(da),'_dx_',num2str(dx),'_','_dm_',num2str(dm),'_','r',num2str(vv),'p',num2str(pp),'_',num2str(S0(1)),'_tr0_',num2str(teta_r0),'_ts0_',num2str(teta_s0),'_a0_',num2str(alpha0),'_x0_',num2str(x0),'_m0_',num2str(m0),'.mat');
fn=fullfile('C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\PM_matfiles',filename);
save(fn)

 

end
