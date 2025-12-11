%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%
% Q. Tao and W. Yuan, "Capacity Analysis of Single-User MIMO-OTFS Systems," in IEEE Transactions on Vehicular Technology, doi: 10.1109/TVT.2025.3588218.
clear all
clc
tic;
P_TR =5;%number of path

M = 10;
N = 10;
N_T = 8;
N_R = 8;
F_M = 1/sqrt(M)*dftmtx(M);
F_N = 1/sqrt(N)*dftmtx(N);
Ptx = eye(M);
Prx = eye(M);
omega = exp(1j*2*pi/(M*N));
Mont=10^1;


%%%%%%  generate channel H   %%%%%%%%%%%%
for j=1:Mont
    
    Delay_TR = randi([0,10],1,P_TR);%integer delay
    Doppler_TR = 2*rand(1,P_TR);%fractional Doppler
    theta_TR_A=1*rand(1,P_TR);
    theta_TR_D=1*rand(1,P_TR);
    h_TR(1) = sqrt(1/2)*(randn(1,1)+ 1j*randn(1,1));%LOS
    h_TR(2:P_TR) = 0.1*sqrt(1/2)*(randn(1,P_TR-1)+ 1j*randn(1,P_TR-1));%NLOS
    H=zeros(M*N*N_R,M*N*N_T);
    F=kron(F_N,eye(M));
    for p = 1:P_TR
        H_p_A(:,:,p)= steeringVec(N_R,theta_TR_A(p))*transpose(steeringVec(N_T,theta_TR_D(p)));
        G(:,:,p)=F*(diag(omega.^((0:M*N-1)*Doppler_TR(p))))...
            *circshift(eye(M*N),Delay_TR(p))*F';
        gg(p)=det(eye(M*N)+G(:,:,p));
        H = H + h_TR(p)*kron(H_p_A(:,:,p),G(:,:,p));
    end
    %%%%%%%%%   SVD  %%%%%%%%%%%%%
    R=rank(H);
    [U,SS,V]=svd(H);
    
    %%%%%%   water filling  %%%%%%%%%
    S=diag(SS);
    lambda=(S(1:R)').^2;
    
    %%%%%% Capacity  %%%%%%%%
    SNR=0:2:20;
    for i=1:length(SNR)
        snr=10^(SNR(i)/10);
        p=waterfill(snr,1./lambda);
        Cap(i,j)=0;
        la=0;
        for r=1:R
            Cap(i,j)=Cap(i,j)+1/(M*N)*log2(1+lambda(r)*p(r));
            la=la+lambda(r)*p(r);
        end
        
        
        %%%%%%%%  Asymptotic Capcity %%%%%%%%%%%%%%%
        Cap_approx(i,j)=log2(1+snr/(M*N)*abs((h_TR(1)))^2*N_T*N_R);
        
    end
    
end


plot(SNR, mean(real(Cap),2),'bv-','Linewidth', 1)
hold on
plot(SNR, mean(real(Cap_approx),2),'K^--','Linewidth', 1)
hold on

h=legend('Capacity, $P=10$','Asymptotic Capcity, $P=10$','Capacity, $P=5$','Asymptotic Capcity, $P=5$','Capacity, $P=1$','Asymptotic Capcity, $P=1$');%,'Capacity LB')
hold on
set(h,'Interpreter','latex')
xlabel('SNR (dB)');
ylabel('Capacity (bit/s/Hz)');
grid on
toc
