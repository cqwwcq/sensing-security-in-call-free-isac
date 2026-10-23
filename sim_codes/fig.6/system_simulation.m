clear all;clc;close all;
xtarget = -75;
ytarget=-xtarget;


PosTarget = [xtarget,ytarget];
if(xtarget<0)
Fake_Target_Pos=PosTarget + [80, -100];
else
Fake_Target_Pos=PosTarget - [80, -100];
end

Ntx = 8; %number of transmit APs
Mtx = 16; %number of antennas at all the APs
Nrx = 4; %number of receive APs
Nue = 2; %number of users
Mue = [16]; %number of antennas at the users

% 优化问题参数
P_T = 50; % dBm

Gamma_sens = 10; % dB  
Gamma_comm = 10; % dB

sigma_noise = -94; % dBm
sigma_ap = -94; % dBm
epsilon = 1e-10;
fprintf('系统配置: %d APs, %d 用户, %d 天线/AP\n', Ntx, Nue, Mtx);
% % 调用CVX求解器
% [X_i, X_t, Z, tau_opt, solve_time] = solve_privacy_isac_bisection_sdp(...
%     H, Txbeamsteering, betark, sigmarcs, ...
%     Ntx, Nrx, Nue, Mtx, P_T, Gamma_sens, Gamma_comm, sigma_noise, sigma_ap);
% 模拟参数

num_runs = 100 % 运行次数（建议 20-50，根据机器性能调整）

% 存储结果的 table

fprintf('开始模拟 (%d 次)...\n', num_runs);
ki=2;

rouMat=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5];


for irou=1:length(rouMat)  
    rou=rouMat(irou);
for itery=1:length(Mue)
HMCx = zeros(Mue(itery),Mtx*Ntx,Nue,num_runs);
miu0=0;
min_sinr=0;
power_total_1=0;
power_total_2=0;
mean_interf_UE=0;
mean_interf_sense=0;
for run = 1:num_runs
    % rng(30);
    fprintf('\n第%i个功率变量，--- 运行 %d/%d ---\n',irou, run, num_runs);
    
    % Step 1: 生成 H 和其他参数
    [lambdac, betark, Dsn, PosAPs, Posue, PosRx, distUe, distTx, distRx, ...
     Txbeamsteering, Rxbeamsteering, Userbeamsteering,Fake_Txbeamsteering,alphark, sigmarcs, H,Z_shd] = ...
     Parameters(Ntx, Nrx, Mtx, Nue, Mue(itery), [xtarget, ytarget],Fake_Target_Pos);
    
 

    
    % Step 3: 运行优化





        fprintf('开始优化...\n');
        if ki==1;

        [gamma_t_opt, w_i_opt, w_t_opt,u_opt, solve_time, status,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense] = run_optimization(...
        Ntx, Nue, Mtx,Nrx,Mue, P_T, Gamma_comm, sigma_noise, sigma_ap,epsilon, ...
        H, Txbeamsteering, Userbeamsteering, sigmarcs,betark,rou);
        elseif ki==2
        [gamma_t_opt, w_i_opt, w_t_opt,u_opt,z_opt, solve_time, status,power_z,miu_z,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense] = AN_run_optimization(...
          Ntx, Nue, Mtx, Nrx,Mue(itery),P_T, Gamma_comm, sigma_noise, sigma_ap,epsilon, ...
          H, Txbeamsteering, Fake_Txbeamsteering, sigmarcs,betark,rou);
        end

    
        if ~status % 如果 status 为 false (0)
            fprintf('优化失败或未找到可行解。\n');
            gamma_t_db = NaN;
            gamma_t_opt = NaN;
            w_i = {}; w_t = {};
        else
            gamma_t_db = 10 * log10(gamma_t_opt);
            fprintf('优化结果: gamma_t=%.6f (%.6f dB) | 时间=%.6f s\n', ...
                    gamma_t_opt, gamma_t_db, solve_time);


  
            end           
                      
            % --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%互信息评估模块结束%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---          
        % end

w_t_opt=reshape(w_t_opt, Mtx, 1, Ntx);
Wx(:,:,:,run) = cat(2, w_i_opt, w_t_opt);
HMCx(:,:,:,run) = H;
Ux=zeros(Mue(itery),length(u_opt),num_runs);
for numu=1:length(u_opt)
    Ux(:,numu,run)=u_opt{numu};
end

if ki==2
miu0=miu0+miu_z;
Zx(:,:,run)= z_opt;
end

min_sinr(run)=min_sinr_dB;
power_total_1=power_total_1+power_1;
power_total_2=power_total_2+power_2;
mean_interf_UE=mean_interf_UE+interf_UE;
mean_interf_sense=mean_interf_sense+interf_sense;
end

if ki==2
miu0=miu0/num_runs;
Zy(:,itery,:,:)=Zx;
end

mean_interf_UE=mean_interf_UE/num_runs;
mean_interf_sense=mean_interf_sense/num_runs;
power_total_1=power_total_1/num_runs;
power_total_2=power_total_2/num_runs;
HMCy((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,:)=HMCx;
Wy(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,:)=Wx;
Uy((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:)=Ux;
end

if ki==2
miu(irou)=miu0;
Z(:,:,:,:,irou)=Zy;
end
HMC(:,:,:,:,irou)=HMCy;
W(:,:,:,:,irou)=Wy;
U(:,:,:,irou)=Uy;
end

if ki==1
save HMC2.mat HMC
save WMC2.mat W
save UMC2.mat U
elseif ki==2
save HMC0.mat HMC
save WMC0.mat W
save ZMC0.mat Z
save UMC0.mat U
end

