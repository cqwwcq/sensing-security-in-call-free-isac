function [lambdac,betark,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,Userbeamsteering,Fake_Txbeamsteering,alphark,sigmarcs,H,Z_shd]=Parameters(Ntx,Nrx,Mtx,Nue,Mue,PosTarget,Fake_Target_Pos)
lambdac = 3e8/1.9e9;
fc=1.9e9;%carrier wavelength
[cn,sn]=signalgeneration(Nue,1); %generating symbols for radar and communication
% Dsn = [sn;cn]; 
Dsn = [cn;sn];%diagonal matrix containing the sensing and communication symbols
%PosAPs = -500 + (500 + 500).*rand(Ntx,2); %Position of transmit APs
%Posue = -500 + (500 + 500).*rand(Nue,2); %Position of Users
% PosAPs = [-500,-500;500,500;-500,500;500,-500;0,-500;0,500;-500,0;500,0]; %location of APs considering 8 Ntx
PosAPs = [-200,-200;200,200;-200,200;200,-200;0,-200;0,200;-200,0;200,0];
%  PosAPs = [
% -437.5 -375; -312.5 -375; -187.5 -375; -62.5 -375; 62.5 -375; 187.5 -375; 312.5 -375; 437.5 -375;
% -437.5 -125; -312.5 -125; -187.5 -125; -62.5 -125; 62.5 -125; 187.5 -125; 312.5 -125; 437.5 -125;
% -437.5  125; -312.5  125; -187.5  125; -62.5  125; 62.5  125; 187.5  125; 312.5  125; 437.5  125;
% -437.5  375; -312.5  375; -187.5  375; -62.5  375; 62.5  375; 187.5  375; 312.5  375; 437.5  375
% ];

%  PosAPs =  [
% -300 -300; -100 -300; 100 -300; 300 -300;
% -300 -100; -100 -100; 100 -100; 300 -100;
% -300  100; -100  100; 100  100; 300  100;
% -300  300; -100  300; 100  300; 300  300
% % ];
%  PosAPs = [
% -100 -100; -33.3 -100; 33.3 -100; 100 -100;
% -100  -33.3; -33.3  -33.3; 33.3  -33.3; 100  -33.3;
% -100   33.3; -33.3   33.3; 33.3   33.3; 100   33.3;
% -100  100; -33.3  100; 33.3  100; 100  100
% ];
%PosAPs = [0,-500;0,500;-500,0;500,0]; %location of APs considering 4 Ntx
%PosAPs = [0,-500;0,500]; %location of APs considering 2 Ntx
% Posue = [300,300;-300,-300;-300,300;300,-300];
% Posue = [100,100;-100,-100;-100,100;100,-100];
Posue = [50,50;-50,-50];
% Posue = [50,50;-50,-50;-50,50;50,-50];

%PosTarget = [25,-350];
%PosRx = -15 + (15 + 15).*rand(Nrx,2);
PosRx = [150,150;-150,150;150,-150;-150,-150]; %location of Receive APs for Nrx= 4
% PosRx = [20,20;-20,20;20,-20;-20,-20]; %location of Receive APs for Nrx= 4
% PosRx = [-250,250;250,-250]; %location of Receive APs for Nrx = 2
%PosRxaux = repmat(PosRx,Ntx-Nrx,1);
%Posueaux = repmat(Posue,Ntx,1);
%PosAPaux = repmat(PosAPs,Nue,1);
%Rrk = corr2(PosRxaux,PosAPs);
%distTR = sqrt((PosRxaux(:,1)-PosAPs(:,1)).^2+(PosRxaux(:,2)-PosAPs(:,2)).^2);

%distUe =  sqrt((Posueaux(:,1)-PosAPs(:,1)).^2+(Posueaux(:,2)-PosAPs(:,2)).^2);
for i=1:Ntx
    for j = 1:Nue
        distUe(i,j) =  sqrt((Posue(j,1)-PosAPs(i,1)).^2+(Posue(j,2)-PosAPs(i,2)).^2);
    end
end

[angleTx, distTx] = compute_angle_dist(PosAPs, PosTarget); %azimuth angle calculation with transmit APs
[angleRx, distRx] = compute_angle_dist(PosRx, PosTarget); %azimuth angle calculation with receive APs
[angleUA, ~] = compute_angle_dist(PosAPs, Posue);
[angle_Fake_Tx,~] = compute_angle_dist(PosAPs,Fake_Target_Pos);
Txbeamsteering = beamsteering(angleTx.', Mtx); % (Ntarget x Ntx x Mtx) %transmit beamsteering
Rxbeamsteering = beamsteering(angleRx.', Mtx); % (Ntarget x Nrx x Mtx) %receive beamsteering
Userbeamsteering = beamsteering(angleUA.', Mtx);% (Nue x Ntx x Mtx) % 
Fake_Txbeamsteering = beamsteering(angle_Fake_Tx.', Mtx);% (Nue x Ntx x Mtx) %
sigmarcs = sqrt(10^(10/10));%单位:平方米

alphark = (sigmarcs).*(randn(1,1)+1i*randn(1,1));

betark = lambdac^2./((4*pi)^3*(distTx.^2)*((distRx.^2).'));
% %Channel coefficient users 
%%%%%%Create the MxK correlated shadowing matrix %%%%%%%
Hb = 15; % Base station height in m
Hm = 1.65; % Mobile height in m
f = 1900; % Frequency in MHz
aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;
phi = 3;
BETAA = zeros(Ntx,Nue);
%%%%M correlated shadowing cofficients of M APs:
d0=0.01;%km
d1=0.05;%km
Dist=zeros(Ntx,Ntx);%distance matrix
Cor=zeros(Ntx,Ntx);%correlation matrix
sigma_shd=8; %in dB
D_cor=0.1;
for m1=1:Ntx
    for m2=1:Ntx
        Dist(m1,m2) = norm(PosAPs(m1,:)-PosAPs(m2,:))*1e-3;
        Cor(m1,m2)=exp(-log(2)*Dist(m1,m2)/D_cor);
    end
end
A1 = chol(Cor,'lower');
x1 = randn(Ntx,1);
sh_AP = A1*x1;
for m=1:Ntx
    sh_AP(m)=(1/sqrt(2))*sigma_shd*sh_AP(m)/norm(A1(m,:));
end

%%%%K correlated shadowing matrix of K terminal:
Dist=zeros(Nue,Nue);%distance matrix
Cor=zeros(Nue,Nue);%correlation matrix

for k1=1:Nue
    for k2=1:Nue
        Dist(k1,k2)=norm(Posue(k1,:)-Posue(k2,:))*1e-3;
        Cor(k1,k2)=exp(-log(2)*Dist(k1,k2)/D_cor);
    end
end
A2 = chol(Cor,'lower');
x2 = randn(Nue,1);
sh_Ter = A2*x2;
for k=1:Nue
    sh_Ter(k)=(1/sqrt(2))*sigma_shd*sh_Ter(k)/norm(A2(k,:));
end

%%% The shadowing matrix:
Z_shd=zeros(Ntx,Nue);
for m=1:Ntx
    for k=1:Nue
        Z_shd(m,k)= sh_AP(m)+ sh_Ter(k);
    end
end

for m=1:Ntx 
    for k=1:Nue
        if distUe(m,k)*1e-3<d0
            betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
        elseif ((distUe(m,k)*1e-3>=d0) && (distUe(m,k)*1e-3<=d1))
            betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(distUe(m,k)*1e-3);
        else
            betadB = -L - 35*log10(distUe(m,k)*1e-3) +Z_shd(m,k);%large-scale in dB
            % betadB = -L - 35*log10(distUe(m,k)*1e-3)+10*randn(1,1);%large-scale in dB

        end
BETAA(m,k)=10^(betadB/10); 
    end
end


phi=3;
for i = 1:Ntx
    % Hax = (sqrt(distUe(i,:).^(-phi)/2))'.*(randn(Nue,Mtx,Mue)+1i*randn(Nue,Mtx,Mue)); %Rayleigh modelling
    % loss = 10.^((36.7+22.7*log10(distUe(i,:)) +26*log10(fc))/10)
    % Hax=1./sqrt(loss').*(randn(Nue,Mtx,Mue)+1i*randn(Nue,Mtx,Mue))./sqrt(2);
    Hax=sqrt(BETAA(i,:)').*(randn(Nue,Mtx,Mue)+1i*randn(Nue,Mtx,Mue))./sqrt(2);
    Husers = permute(Hax,[3 2 1]);
    H(:,(1+(i-1)*Mtx):i*Mtx,:) = Husers; %Number of receive antennas X Number of transmit antennas x Number of BS X Number of users
end
end 