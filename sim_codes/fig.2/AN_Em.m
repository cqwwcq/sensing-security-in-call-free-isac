%clear all
%clc
function [Rall]= AN_Em(Dsn,nvariance,Ntx,Nue,Mue,H,sigma2,W,Z,Mtx)
Covnd = nvariance;
cn = Dsn(1:Nue,:);
sn = Dsn(Nue+1,:);
for temp = 1:128
jam_symbol(temp) = exp(1j * 2 * pi * rand()); 
end
for k = 1:Ntx
for temp = 1:128
        yt = H(:,(1+(k-1)*Mtx):k*Mtx,1)*W(:,Nue+1,k)*sn(:,temp);%size of W is[Mtx Nue+1 Ntx]
        Waux = W;
        Waux(:,[1,Nue+1],:) = [];
        % Waux(:,[Nue+1],:) = [];
        cnaux = cn;
        cnaux(1,:) = [];
        yue =  H(:,(1+(k-1)*Mtx):k*Mtx,1)*Waux(:,:,k)*cnaux +H(:,(1+(k-1)*Mtx):k*Mtx,1)*Z(:,k)*jam_symbol(temp)+ sigma2;       
        y = yue(:,temp)+yt;
        sdhataux = pinv(H(:,(1+(k-1)*Mtx):k*Mtx,1)'*H(:,(1+(k-1)*Mtx):k*Mtx,1))*H(:,(1+(k-1)*Mtx):k*Mtx,1)'*y;
        cvg = norm(W(:,Nue+1,k)*sn(:,temp) + Waux(:,:,k)*cnaux);
        Eh = H(:,(1+(k-1)*Mtx):k*Mtx,1)+ (normrnd(0,1e-2,Mue,Mtx)+1i.*normrnd(0,1e-2,Mue,Mtx));     
        Omegah=inv(sdhataux*sdhataux'+eye(Mtx)/(1e-2));
 %EM algorithm
    while cvg >= 1e-5
        a = (10 + 2*Mue)/2;
        Cl = ((y - Eh*(sdhataux))'*(y- Eh*(sdhataux))+sdhataux'*Omegah*sdhataux*Mue)/Covnd;
        b = (10 +2* Cl)/2;
        qu = gamrnd(a,b);
        Eu = a./b;
        Omegah = inv(eye(Mtx)/(1e-2)+ (sdhataux)*(sdhataux)'.*(real(Eu)/Covnd));
        Eh = (conj(Omegah)*(H(:,(1+(k-1)*Mtx):k*Mtx,1)'/(1e-2) + (sdhataux)*y'*(Eu/Covnd)))';
        %qh = normrnd(0,sqrt(potNoiseEst),Mt,1,N);
        Vh = chol(Omegah*Mue);
        Htil = [Eh;Vh];
        ytil = [y;zeros(Mtx,1)];
       %funs = @(Wrstil)norm(ytil - Htil*Wrstil)^2;
       %mini = [];
       %maxi = [];
       %sdhat = fmincon(funs,sdhataux,mini,maxi);
       sdhat = inv(Htil'*Htil)*Htil'*ytil;
       cvg = norm(sdhat-sdhataux);
       sdhataux = sdhat;
    end
    Raux(:,:,temp) = sdhataux*sdhataux';
end
    Rall(:,:,:,k) = Raux; 
end
% for k = 1:Ntx
%      thetal = linspace(angsAPs(k,1),angsAPs(k,2),901); %angle grid sample
%      figure
%      plot(rad2deg(thetal),pt(:,:,k))
%  end

end
