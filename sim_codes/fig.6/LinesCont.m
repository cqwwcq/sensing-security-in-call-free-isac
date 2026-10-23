 clear all
clc
PosTarget = [-75,75]; %true position of the target
Fake_Target_Pos = PosTarget + [80, -100]; 
P_T = 50; % dBm
P_T_lin = 10^(P_T/ 10) / 1000;      % W
% ratio_AN=0.5984;
% [Ntx,Nrx,Mtx,Nue,gammac,nvariance,lambdac,Dsn,PosAPs,Posue,PosRx,distTx,distRx,Txbeamsteering,Rxbeamsteering,alphark,sigmarcs,H,sigma2,W]=ParametersCont(PosTarget);
% rng(53);
 % rng(56);
% rng(25);
% rng(21);
% rng(14);
% rng(35);
% rng(36)
% rng(39);
% rng(37);
% rng(47);
Ntx = 8; %number of transmit APs
Nrx = 4; %number of receive APs
Mtx = 16; %number of antennas at all the APs
Nue = 2; %number of users
Mue = [16];
numrun=100;
itery = 1;
 % gammac = 10^(3/10); %threshold for the communication
gammac = 10^(3/10)/1000;
 % nvariance = 10^(-94/100)/1000;
nvariance = 10^(-94/10)/1000;
% sigma2 = abs(nvariance.*randn(1,1));
% seed=[53 56 25 21 13 35 36 39 37 47]
seed= randi([0, 2^32-1],1,10);
tic


rouMat=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5];

ki=2;
for irou=1:length(rouMat)
    rou=rouMat(irou);

    if ki==1
        HMC=load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\HMC2');
        W  =load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\WMC2');
        U  =load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\UMC2');
    elseif ki==2       

        HMC=load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\HMC0');
        W  =load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\WMC0');
        Z  =load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\ZMC0');
        U  =load('C:\Users\pc\Desktop\代码\AN_simulation\fig.6\UMC0');

    end
HMCy=HMC.HMC;
Wy=W.W;
Uy=U.U;
if ki==2 
    Zy=Z.Z;
    Z=Zy(:,:,:,:,irou);
end
HMC=HMCy(:,:,:,:,irou);
W=Wy(:,:,:,:,irou);
U=Uy(:,:,:,irou);





xcont = -195:20:195; %center of the search area cells
ycont = 195:-20:-195;
% xcont = -475:50:475; %center of the search area cells
% ycont = 475:-50:-475;
%n = length(xcont) ;
%m = length(ycont) ; 
[Xcont,Ycont] = meshgrid(xcont,ycont) ; 
x =-205:20:205; %Entire search area
y = 205:-20:-205;
% x =-500:50:500; %Entire search area
% y =500:-50:-500;
n = length(x) ;
m = length(y) ; 
[X,Y] = meshgrid(x,y) ; 

 %modify this part indication the number of antennas at the UE or modify accordingly to evaluate other parameter.
if Ntx ==8 %interval of the steering vectors 
% angsAPs = [0,pi/2;-pi,-pi/2;-pi/2,0;pi/2,pi;-pi/2,pi/2;-pi/2,pi/2;pi/2,-pi/2;-pi/2,pi/2];
angsAPs = [0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;];
elseif Ntx == 4
angsAPs = [-pi/2,pi/2;-pi/2,pi/2;-pi/2,pi/2;-pi/2,pi/2];
elseif Ntx == 16
angsAPs = [0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi];

else
angsAPs = [pi/2,-pi/2;pi/2,-pi/2];
end

%for itery = 1:length(Mue)
for s=1:length(seed)
pd = 0;
pa = 0;
Poscont = zeros(100,2);
rng(seed(s));
sigma2 = sqrt(nvariance/2) * (randn(Mue(itery),1) + 1j*randn(Mue(itery),1));
[lambdac,betark,Dsn,PosAPs,Posue,PosRx,distUe,distTx,distRx,Txbeamsteering,Rxbeamsteering,Userbeamsteering,Fake_Txbeamsteering,alphark,sigmarcs,H,Z_shd]=Parameters(Ntx,Nrx,Mtx,Nue,Mue(itery),PosTarget,Fake_Target_Pos);

for MC = 1:numrun %obtaining the covariance matrix with the EM algorithm
    if ki==1
    Z_AN = generate_projected_AN(Nue,Mtx,Ntx,HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC), U((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,MC), P_T_lin, rou);
    Rall(:,:,:,:,MC)= AN_Em(Dsn,nvariance,Ntx,Nue,Mue(itery),HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC),sigma2,W(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC),Z_AN,Mtx);
    % Rall(:,:,:,:,MC)= Em(Dsn,nvariance,Ntx,Nue,Mue(itery),HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC),sigma2,W(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC),Mtx);
    elseif ki==2
    % Z_AN = generate_fakeangle_projected_Z(Nue,Mtx,Ntx,HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC), U((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,MC), P_T_lin, ratio_AN,PosAPs,Fake_Target_Pos);
    Z_AN=Z(:,itery,:,MC);
    Rall(:,:,:,:,MC)= AN_Em(Dsn,nvariance,Ntx,Nue,Mue(itery),HMC((1+(itery-1)*(Mue(itery))):itery*(Mue(itery)),:,:,MC),sigma2,W(:,(1+(itery-1)*(Nue+1)):itery*(Nue+1),:,MC),Z_AN,Mtx);
    end
end
    %Rtotal(:,:,:,:,(1+(itery-1)*100):itery*100) = Rall;
%end
%figure
for MC = 1:numrun %verification for each MC iteration of the probability of detection
angmax = zeros(1,Ntx);


ptotal = zeros(3604,128);
pt = zeros(1,3604,Ntx);
for k =1:Ntx
    %if abs(PosAPs(k,1)) == abs(PosAPs(k,2))
        thetal = linspace(angsAPs(k,1),angsAPs(k,2),3604); %angle grid sample
        a = (1i*2*pi*(1/2)*(1:1:(Mtx-1)).*sin(thetal'))'; %beam steering
        %atargets = (i*2*pi*(1/2)*(1:1:(M-1)).*sin([theta2 theta1 theta3]'))';
        asteer = [ones(1,length(thetal));exp(a)];
        py = zeros(1,length(thetal));
        % for temp = 1:128
        % for i=1:length(thetal)
        %     py(i) = real(asteer(:,i)'*Rall(:,:,temp,k,MC)*asteer(:,i)); %beampattern
        % end
        % ptotal(:,temp) = py;    
        % end    
        % pt(:,:,k) = sum(ptotal,2)'; %beampattern for all the APs
        R_sum = sum(Rall(:,:,:,k,MC), 3);   % Mtx x Mtx
        % 批量计算：a_i' * R_sum * a_i 对所有i
        % = diag(asteer' * R_sum * asteer)
        temp_mat = asteer' * R_sum;           % 3604 x Mtx
        py_vec   = real(sum(temp_mat .* asteer.', 2))';  % 1 x 3604
        pt(:,:,k) = py_vec;
end
contgrid = zeros(length(x)-1,length(y)-1,Ntx);
[maxim,ind] = max(pt); %get the index and value of the maximum point of the beampattterns
% xreta = -200:0.0004:200; 
xreta = -200:0.001:200;
% xreta = -500:0.001:500; 
yreta = zeros(1,length(xreta),Ntx);
% angcont = [0,pi/2;-pi,-pi/2;-pi/2,0;pi/2,pi;0,pi;0,-pi;pi/2,-pi/2;-pi/2,pi/2];
angcont = [0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;];

% angcont = [0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi;0,2*pi];

for k = 1:Ntx%???????????????????
    thetal = linspace(angcont(k,1),angcont(k,2),3604); %angle grid sample
    angmax(k) = thetal(ind(:,:,k));
    yreta(:,:,k) = tan(angmax(k)).*(xreta-PosAPs(k,1))+PosAPs(k,2); %calculate the direction line
for i = 1:m-1
  for j = 1:n-1
    % Get the points of each grid 
    P = [X(i,j) Y(i,j) ; 
         X(i,j+1) Y(i,j+1) ; 
         X(i+1,j+1) Y(i+1,j+1) ;
         X(i+1,j) Y(i+1,j)] ; 
         
    idx = any(inpolygon(xreta,yreta(:,:,k),P(:,1),P(:,2)));
    contgrid(i,j,k) = idx;
    %iwant = [x(idx) y(idx) z(idx)] ; 
    %plot(x(idx),y(idx),'.')
    %drawnow
   end
end
end

contTotal = (sum(contgrid,3)); %obtain a matrix with the count of cells that are crossed by the lines
maximum = max(max(contTotal)); %verify the line with the maximum number of lines crossed.
[xc,yc]=find(contTotal==maximum);
if length(xc)>1
    aux = randi(length(xc)); 
    xf = xc(aux);
    yf = yc(aux);
    Poscont(MC,:) = [Xcont(xf,yf),Ycont(xf,yf)];
else
    xf = xc;
    yf = yc;
    Poscont(MC,:) = [Xcont(xf,yf),Ycont(xf,yf)];
end
if Poscont(MC,:) == PosTarget %calculate the probability of detection?????????????
    pd = pd + 1;
end
end
fprintf("第%i个功率变量，第%i次实验，检测目标%i个\n",irou,s,pd);
final_Poscont(:,:,s,irou)=Poscont;
final_pd(s,irou)=pd;
end
finalans(irou)=sum(final_pd(:,irou))/length(final_pd(:,irou));


end
% Posconttotal = [Posconttotal Poscont];
toc
% save Posconts1knownAPsm81.mat Poscont
