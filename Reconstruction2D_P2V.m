%-------------------------------------------------------------------------
%    This code reads the distribution in position space
%        and reconstructs the fre-energy landscape
%-------------------------------------------------------------------------

% ---- k-space parameters ----
M=19;MM=1001; [n,m]=deal(-floor(M/2):floor(M/2));
% ---- potential parameters
kBTx=1; kBTy=1; gx=1; gy=1; L=2;

%----------------  Pkl from P_pos 

Pkl=fftshift(fft2(real(P_pos))/(M^2)); % Fourier Series 
Inxnorm=floor(size(P_pos)/2)+1;          % For Normalization

%------------ This get rid of unused frequencies 

Pkl=Pkl/Pkl(Inxnorm(1),Inxnorm(2));      
Pkl=Pkl(Inxnorm(1)-floor(M/2):Inxnorm(1)+floor(M/2),...
    Inxnorm(2)-floor(M/2):Inxnorm(2)+floor(M/2)); 
    
%------------ This reconstructs the potential

[KK,BB,P]=InversionMatrix(Pkl,M,fx,fy,L);
[VV]=SolverP2V(BB,KK,P,M);
[Vrec,~,~]=kspace2position1(VV,n,MM,0,L);
Vrec=real(Vrec)-max(real(Vrec(:))); %fix reference level
