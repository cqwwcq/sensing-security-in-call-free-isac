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

% Optimization problem parameters
P_T = 50; % dBm

Gamma_sens = 10; % dB
Gamma_comm = 10; % dB

sigma_noise = -94; % dBm
sigma_ap = -94; % dBm
epsilon = 1e-10;
fprintf('System configuration: %d APs, %d users, %d antennas/AP\n', Ntx, Nue, Mtx);
% % Call the CVX solver
% [X_i, X_t, Z, tau_opt, solve_time] = solve_privacy_isac_bisection_sdp(...
%     H, Txbeamsteering, betark, sigmarcs, ...
%     Ntx, Nrx, Nue, Mtx, P_T, Gamma_sens, Gamma_comm, sigma_noise, sigma_ap);
% Simulation parameters
num_runs = 1 % Number of runs (recommended: 20-50, depending on machine performance)

% Table for storing results
fprintf('Starting simulation (%d runs)...\n', num_runs);
ki=2;

if Mtx==16
    rouMat=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5];
    % rouMat=[0.1039, 0.1879, 0.2565, 0.3163, 0.3665, 0.4065, 0.4475, 0.4795, 0.5072, 0.5351, 0.5582, 0.5827...
    %     , 0.5984, 0.6144, 0.6308, 0.7, 0.75, 0.8, 0.85, 0.88];
elseif Mtx==8
    rouMat=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5];
    % rouMat=[0.0019, 0.1235, 0.2234, 0.3049, 0.3667, 0.4129, 0.4612,0.4976,...
    % 0.5322, 0.5588, 0.5811, 0.6064, 0.6252, 0.6444, 0.6644, 0.6748, 0.7, 0.75, 0.8, 0.85, 0.9];
end

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
            fprintf('\nPower variable %i, --- run %d/%d ---\n',irou, run, num_runs);
            % Step 1: Generate H and other parameters
            [lambdac, betark, Dsn, PosAPs, Posue, PosRx, distUe, distTx, distRx, ...
                Txbeamsteering, Rxbeamsteering, Userbeamsteering,Fake_Txbeamsteering,alphark, sigmarcs, H,Z_shd] = ...
                Parameters(Ntx, Nrx, Mtx, Nue, Mue(itery), [xtarget, ytarget],Fake_Target_Pos);
            % Step 3: Run optimization
            fprintf('Starting optimization...\n');
            if ki==1;
                [gamma_t_opt, w_i_opt, w_t_opt,u_opt, solve_time, status,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense] = run_optimization(...
                    Ntx, Nue, Mtx,Nrx,Mue, P_T, Gamma_comm, sigma_noise, sigma_ap,epsilon, ...
                    H, Txbeamsteering, Userbeamsteering, sigmarcs,betark,rou);
            elseif ki==2
                [gamma_t_opt, w_i_opt, w_t_opt,u_opt,z_opt, solve_time, status,power_z,miu_z,power_total,min_sinr_dB,power_1,power_2,interf_UE,interf_sense] = AN_run_optimization(...
                    Ntx, Nue, Mtx, Nrx,Mue(itery),P_T, Gamma_comm, sigma_noise, sigma_ap,epsilon, ...
                    H, Txbeamsteering, Fake_Txbeamsteering, sigmarcs,betark,rou);
            end
            if ~status % If status is false (0)
                fprintf('Optimization failed or no feasible solution was found.\n');
                gamma_t_db = NaN;
                gamma_t_opt = NaN;
                w_i = {}; w_t = {};
            else
                gamma_t_db = 10 * log10(gamma_t_opt);
                fprintf('Optimization result: gamma_t=%.6f (%.6f dB) | time=%.6f s\n', ...
                    gamma_t_opt, gamma_t_db, solve_time);
            end
            % --- Mutual information evaluation module ends ---
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
