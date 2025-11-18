% function [t,q_bot]=reading_flx_Hy
function [t,vtop,q_bot,inf,Ev]=reading_flx_Hy
% function [t,vtop,q_bot,inf,Ev,t_m,P]=reading_flx_Hy

work_dir='C:\Hydrus_test\simulation\';
% work_dir='C:\Hydrus_test\sim_par_err\';
% mkdir(work_dir);

num_sim=1;
path=cell(1,num_sim);
k=1;
% for k=1:num_sim
path{k}=strcat(work_dir,'sim_param_',num2str(k));
% path{k}=strcat(work_dir,'sim_param_err_',num2str(k));
fname=fullfile(path{k},'T_level.out');
% reading output
fluxes=importdata(fname);
t=fluxes.data(:,1);
vtop=fluxes.data(:,4);
q_bot=fluxes.data(:,6);
inf=fluxes.data(:,18);
Ev=fluxes.data(:,19);

% fname2=fullfile(path{k},'Meteo.out');
% % reading output
% meteo=importdata(fname2);
% t_m=meteo.data(:,1);
% P=meteo.data(:,9)/10;


end

