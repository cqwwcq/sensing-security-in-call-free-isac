function [Rall, cap_hits] = AN_Em(Dsn,nvariance,Ntx,Nue,Mue,H,sigma2,W,Z,Mtx)
% Adversary-side EM reconstruction in the presence of AN, as in Section II-D.
% Identical to the reference implementation except for an iteration cap on
% the EM loop, which guards the parallel sweep against a non-converging
% realization; cap_hits reports how often it was reached.

MAX_EM_ITER = 200;
cap_hits = 0;
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
        cnaux = cn;
        cnaux(1,:) = [];
        yue =  H(:,(1+(k-1)*Mtx):k*Mtx,1)*Waux(:,:,k)*cnaux +H(:,(1+(k-1)*Mtx):k*Mtx,1)*Z(:,k)*jam_symbol(temp)+ sigma2;
        y = yue(:,temp)+yt;
        sdhataux = pinv(H(:,(1+(k-1)*Mtx):k*Mtx,1)'*H(:,(1+(k-1)*Mtx):k*Mtx,1))*H(:,(1+(k-1)*Mtx):k*Mtx,1)'*y;
        cvg = norm(W(:,Nue+1,k)*sn(:,temp) + Waux(:,:,k)*cnaux);
        Eh = H(:,(1+(k-1)*Mtx):k*Mtx,1)+ (normrnd(0,1e-2,Mue,Mtx)+1i.*normrnd(0,1e-2,Mue,Mtx));
        Omegah=inv(sdhataux*sdhataux'+eye(Mtx)/(1e-2));
 %EM algorithm
    it_em = 0;
    while cvg >= 1e-5 && it_em < MAX_EM_ITER
        it_em = it_em + 1;
        a = (10 + 2*Mue)/2;
        Cl = ((y - Eh*(sdhataux))'*(y- Eh*(sdhataux))+sdhataux'*Omegah*sdhataux*Mue)/Covnd;
        b = (10 +2* Cl)/2;
        qu = gamrnd(a,b); %#ok<NASGU>
        Eu = a./b;
        Omegah = inv(eye(Mtx)/(1e-2)+ (sdhataux)*(sdhataux)'.*(real(Eu)/Covnd));
        Eh = (conj(Omegah)*(H(:,(1+(k-1)*Mtx):k*Mtx,1)'/(1e-2) + (sdhataux)*y'*(Eu/Covnd)))';
        Vh = chol(Omegah*Mue);
        Htil = [Eh;Vh];
        ytil = [y;zeros(Mtx,1)];
       sdhat = inv(Htil'*Htil)*Htil'*ytil;
       cvg = norm(sdhat-sdhataux);
       sdhataux = sdhat;
    end
    if it_em >= MAX_EM_ITER
        cap_hits = cap_hits + 1;
    end
    Raux(:,:,temp) = sdhataux*sdhataux';
end
    Rall(:,:,:,k) = Raux;
end

end
