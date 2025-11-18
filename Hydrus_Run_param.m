% function [t,S]=Hydrus_Run_param(S0,alpha0,K0)
% function [t,S]=Hydrus_Run_param(S0,alpha0,n0,K0)
% function [t,S]=Hydrus_Run_param(S0,alpha0,n0,m0)
% function [t,S]=Hydrus_Run_param(S0,alpha0,x0,m0)
function [t,S]=Hydrus_Run_param(S0,teta_r0,teta_s0,alpha0,x0,m0)
% this function gets S0 and shp as input and its output is the SM profile at
% each timestep.
delete('C:\Hydrus_test\simulation\sim_param_1\*');

z=[0:0.01:0.5,0.52:0.02:1,1.03:0.03:1.51];
lz=length(z);
num_node=lz;

%%% I added this for feasibility analysis (changing # of obs) for previous
%%% results these line did not use. the values based on the range of soil
%%% texture class
if x0<log(1)
    x0=log(1);
end
if x0>log(2.3)
    x0=log(2.3);
end

if alpha0<0
    alpha0=0.001;
end

% if teta_r0<0
%     teta_r0=0.0001;
% end
% if teta_r0>0.1
%     teta_r0=0.1;
% end

if teta_s0<0.2
    teta_s0=0.2;
end
if teta_s0>0.7
    teta_s0=0.7;
end


K0=exp(m0);
n0=exp(abs(x0));
%[teta_r, teta_s, alpha,n,Ks,tau,thm,tha,thk,kk]
Par0=[teta_r0 teta_s0 alpha0 n0 K0 0.5 0.39 0.1 0.39 31.44];  

%%% generated these .dat and .in file by running Hydrus IF, then used the
%%% format and modified with my own values.
hydrus_exec='C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\H1D_calc.exe';
profileDAT='C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\jan4,2021\profile.dat';
selectorIN='C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\jan4,2021\selector.in';
meteoIN='C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\jan4,2021\METEO.IN';
atmosphIN='C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run\jan4,2021\ATMOSPH.in';
work_dir='C:\Hydrus_test\simulation\';
% mkdir(work_dir);

num_sim=1;
path=cell(1,num_sim);
k=1;
% for k=1:num_sim
path{k}=strcat(work_dir,'sim_param_',num2str(k));
% mkdir(path{k});
%copy selector.in from reference directory to the simulation directory
% copyfile(selectorIN,path{k});
%copy METEO.in from reference directory to the simulation directory
copyfile(meteoIN,path{k});
%copy atmosph.in from reference directory to the simulation directory
copyfile(atmosphIN,path{k});
%copy profile.dat from reference directory to the simulation directory
copyfile(profileDAT,path{k});

%%%manipulate profile.dat

filename=fullfile(path{k},'profile.dat');
fileID = fopen(filename);
C = textscan(fileID,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s');
fclose(fileID);

for i=5:97
    C{1,3}{i,1}= num2str(S0(i-4));
end

A = cellfun(@(x) x(1:105),C,'UniformOutput',false);

while any(cellfun(@iscell,A))
    A = [A{cellfun(@iscell,A)} A(~cellfun(@iscell,A))];
end

writecell(A,filename,'Delimiter','space');


%%% manipulate Selector.in for each run
fileID_out=fopen(strcat(path{k},'\SELECTOR.in'),'wt');
fileID_in=fopen(selectorIN);
skip_lines=26;
for l=1:(skip_lines)
    x=fgetl(fileID_in);
    fprintf(fileID_out,'%s\n',x);
end

out_par=Par0;  %out_par(3)=0.005;
fprintf(fileID_out,'%f %f %f %f %f %f %f %f %f %f\n',out_par');
fgetl(fileID_in);
skip_lines_end=11;
for l=1:(skip_lines_end)
    x=fgetl(fileID_in);
    fprintf(fileID_out,'%s\n',x);
end
fclose('all');


% to create level_0.dir
fname=fullfile('C:\Users\parisaheidary\Box\Mywork_Recharge\Hydrus\Hydrus_Run','level_01.dir');
fid = fopen(fname, 'w');
level0DIR=fullfile(path{k});  %%path for run folder which all inputs are there
fprintf(fid, '%s', level0DIR);
fclose(fid);

% run hydrus
system(hydrus_exec);


% reading output
filename1=fullfile(path{k},'Obs_Node.out');
fileID1 = fopen(filename1);
CC = textscan(fileID1,'%s','delimiter','\n');
lcc=length(CC{1});
lper=lcc-11-1; %11 is the skip_line number, 1 is for end
fclose(fileID1);

filename11=fullfile(path{k},'Obs_Node.out');
fid_out = fopen(filename11);
C = textscan(fid_out,'%[^\n]','HeaderLines',11);
fclose(fid_out);
C = [C{:}];
for i=1:size(C,1)-1
    % strfind(C{1,2},'*');
    C{i,1} = strrep(C{i,1},'*','0');
    
end
fname1=fullfile(path{k},'temporary.txt');
filePh = fopen(fname1,'w');
fprintf(filePh,'%s\n',C{:});
fclose(filePh);
fname2=fullfile(path{k},'temporary.txt');
fidd=fopen(fname2);

% reading outputs
%     fileID_out=fopen(strcat(path{k},'\Obs_Node.out'));
%     skip_lines=11;
%     for kk=1:(skip_lines)
%         x=fgetl(fileID_out);
%     end

temp=fscanf(fidd,'%f',[3*num_node+1,lper]);
t=temp(1,:);
S=temp(3:3:3*num_node,:);
fclose(fidd);





end