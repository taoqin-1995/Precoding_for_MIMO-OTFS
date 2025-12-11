%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%
% Q. Tao and W. Yuan, "Capacity Analysis of Single-User MIMO-OTFS Systems," in IEEE Transactions on Vehicular Technology, doi: 10.1109/TVT.2025.3588218.

clear all
tic
clc
P_TR = 10;%number of path
M = 20; 
N = 10;
N_T = 8;
N_R = 8;
F_M = 1/sqrt(M)*dftmtx(M);
F_N = 1/sqrt(N)*dftmtx(N);
Ptx = eye(M);
Prx = eye(M);
omega = exp(1j*2*pi/(M*N));
%%%%%%generate channel H%%%%%%%%%%%%
Mont=10^0;
SNR=0:2:20;
for j=1:Mont
Delay_TR = randi([0,10],1,P_TR);%integer delay 
Doppler_TR = 2*rand(1,P_TR);%fractional Doppler
theta_TR_A=1*rand(1,P_TR);
theta_TR_D=1*rand(1,P_TR);
h_TR(1) = sqrt(1/2)*(randn(1,1)+ 1j*randn(1,1));
h_TR(2:P_TR) = 0.1*sqrt(1/2)*(randn(1,P_TR-1)+ 1j*randn(1,P_TR-1));
H=zeros(M*N*N_R,M*N*N_T);
F=kron(F_N,eye(M));
for p = 1:P_TR
     H_p_A(:,:,p)= steeringVec(N_R,theta_TR_A(p))*(steeringVec(N_T,theta_TR_D(p)))';
     G(:,:,p)=F*(diag(omega.^((0:M*N-1)*Doppler_TR(p))))*...
         circshift(eye(M*N),Delay_TR(p))*F';
     H= H+ h_TR(p)*kron(H_p_A(:,:,p),G(:,:,p));
end

R=rank(H);
min(N_T,min(N_R,P_TR))*M*N;
[U,SS,V]=svd(H); 
%%%%%%%%%water filling%%%%%%%%%%%%%%%
S=diag(SS);
lambda=(S(1:R)').^2;

for i=1:length(SNR)
snr=10^(SNR(i)/10);

%%%%%%%%%  Exact Capacity (Theoretical) %%%%%%%%%
p=waterfill(snr,1./lambda);
Cap(i,j)=0;
for r=1:R
    Cap(i,j)=Cap(i,j)+1/(M*N)*log2(1+lambda(r)*p(r));
end

%%%%%%%%%Exact Capacity (Using optimal beamforming and combining vectors)%%%%%%%%%%%
 W_opt=V(:,1:R)*(diag(p))^(1/2);
 F_opt=U(:,1:R)';
 Cap_accurate(i,j)=1/(M*N)*log2(det(eye(R)+W_opt'*H'*F_opt'*F_opt*H*W_opt));
 
%%%%%%%  LoS at T\&R  %%%%%%%
W_full_TR=sqrt(snr/(M*N*N_R))*kron(steeringVec(N_T,theta_TR_D(1)),...
    ((circshift(eye(M*N),Delay_TR(1)))*F')');
F_full_TR=sqrt(1/(N_T))*kron((steeringVec(N_R,theta_TR_A(1)))',...
    (F*diag(omega.^((0:M*N-1)*Doppler_TR(1))))');
Cap_full_TR(i,j)=1/(M*N)*log2(det(eye(M*N)+W_full_TR'*H'*F_full_TR'*F_full_TR*H*W_full_TR));

%%%%%%%%%LoS at R%%%%%%%
W_par_T=sqrt(snr/(M*N*N_R))*kron(steeringVec(N_T,theta_TR_D(1)),...
    (F));
F_par_T=sqrt(1/(N_T))*kron((steeringVec(N_R,theta_TR_A(1)))',...
    (F*(diag(omega.^((0:M*N-1)*Doppler_TR(1)))*circshift(eye(M*N),Delay_TR(1))))');
Cap_par_T(i,j)=1/(M*N)*log2(det(eye(M*N)+W_par_T'*H'*F_par_T'*F_par_T*H*W_par_T));

%%%%%%%%% LoS at T%%%%%%%
W_par_R=sqrt(snr/(M*N*N_R))*kron(steeringVec(N_T,theta_TR_D(1)),...
    (diag(omega.^((0:M*N-1)*Doppler_TR(1)))*(circshift(eye(M*N),Delay_TR(1)))*(F'))');
F_par_R=sqrt(1/(N_T))*kron((steeringVec(N_R,theta_TR_A(1)))',...
    (F)');
Cap_par_R(i,j)=1/(M*N)*log2(det(eye(M*N)+W_par_R'*H'*F_par_R'*F_par_R*H*W_par_R));

%%%%%%  Asymptotic Capcity %%%%%%%%%%
Cap_1(i,j)=log2(1+snr/(M*N)*abs((h_TR(1)))^2*N_T*N_R);

end
end



plot(SNR, mean(real(Cap_accurate),2),'r-','Linewidth', 1) 
hold on
plot(SNR, mean(real(Cap_full_TR),2),'b^-','Linewidth', 1)
hold on
plot(SNR, mean(real(Cap_par_T),2),'gv-','Linewidth', 1)
hold on
plot(SNR, mean(real(Cap_par_R),2),'m*-','Linewidth', 1)
hold on
plot(SNR, mean(real(Cap),2),'y--','Linewidth', 1.5) 
hold on
plot(SNR, mean(real(Cap_1),2),'k--','Linewidth', 1)
hold on
h=legend('Full CSI at T\&R','LoS at T\&R','LoS at R','LoS at T','Exact Capacity','Asymptotic Capacity');%,'Capacity LB')
hold on
set(h,'Interpreter','latex') 
xlabel('SNR (dB)');
ylabel('Capacity (bit/s/Hz)');
grid on
toc
