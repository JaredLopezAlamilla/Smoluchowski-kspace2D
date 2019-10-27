function [Fs]=iterativo2D(fs,P_pos,Vorg,Jxnm,Jynm,L,Ns,MSEs) % 2D
%-------------------------------------------------------------------------
%        RECONSTRUCT  POTENTIAL  FROM  DISTRIBUTION    USING
%                INVERSE MATRIX FOR F-P EQ. 2D
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%                            Initial Parameters
%-------------------------------------------------------------------------
fx=fs(1,1); fy=fs(2,1); % loads force guess 2D
kBTx=1; gx=1; kBTy=1; gy=1; 
M=19; MM=1001; [n,~]=deal(-floor(M/2):floor(M/2));
k=floor(M/2); ks=-k:1:k; % size 2k+1
%-------------------------------------------------------------------------
%    This part of the code reads the distribution in position space
%                 and calculates its k-space form 
%-------------------------------------------------------------------------
Pkl=fftshift(fft2(real(P_pos))/(M^2));   % Fourier Series 
Inxnorm=floor(size(P_pos)/2)+1;          % For Normalization
%------------ This get rid of nonused frequencies 
Pkl=Pkl/Pkl(Inxnorm(1),Inxnorm(2));      
Pkl=Pkl(Inxnorm(1)-floor(M/2):Inxnorm(1)+floor(M/2),...
    Inxnorm(2)-floor(M/2):Inxnorm(2)+floor(M/2)); 
%-------------------------------------------------------------------------
%            BUILDS    SPARE DIAGONAL MATRIX FOR F-P EQ. 2D
%-------------------------------------------------------------------------
W=spdiags((ks.*(ks-0))',0,length(ks),length(ks));
D=spdiags(ones(1,length(ks))',0,length(ks),length(ks));

for l=1:k           
    Dg1=spdiags(ones(1,length(ks))',l,length(ks),length(ks));   
    Dg2=spdiags(ones(1,length(ks))',-l,length(ks),length(ks)); 
    Wg1=spdiags(circshift((ks.*(ks+l)),[0,l])',l,length(ks),length(ks));
    Wg2=spdiags(circshift((ks.*(ks-l)),[0,-l])',-l,length(ks),length(ks));
    
    D=D+Dg1+Dg2;  % This matrix contains the pattern of the full matrix     
    W=W+Wg1+Wg2;
end

D=repmat(D,length(ks));   % size (2k+1)^2 X (2k+1)^2 
Wy=repmat(W,length(ks));  
Wx=[];                    

for l=1:length(ks)
    for r=1:length(ks)    
    Wx(1+(l-1)*length(ks):l*length(ks),1+(r-1)*length(ks):r*length(ks))=W(l,r)*D(1+(l-1)*length(ks):l*length(ks),1+(r-1)*length(ks):r*length(ks));
    end    
end

P=reshape(Pkl,length(ks)^2,1).';
%--------- D2 matrix contains the matrix with P values properly sorted

D2=spdiags(ones(1,length(ks)^2)',0,length(ks)^2,length(ks)^2);

for l=1:floor(length(ks)^2/2)      
    D2g1=spdiags(P(ceil(length(ks)^2/2)-l)*ones(1,length(ks)^2)',l,length(ks)^2,length(ks)^2);    
    D2g2=spdiags(P(ceil(length(ks)^2/2)+l)*ones(1,length(ks)^2)',-l,length(ks)^2,length(ks)^2);     
    D2=D2+D2g1+D2g2;
end
%-------------------------------------------------------------------------
%         Jx and Jy contain the terms k(k-q)Pq which multiplies Vk
%-------------------------------------------------------------------------
Jx=Wx.*D2;
Jy=Wy.*D2;
%-------------------------------------------------------------------------
%         Dx and Dy contain the terms Kx^2+Ky^2 which multiplies Pk
%-------------------------------------------------------------------------
Dy=repmat(ks,1,length(ks));
Dx=sort(repmat(ks,1,length(ks)));

K=spdiags((Dx.*Dx-1j*fx*L*Dx./(2*pi)+Dy.*Dy-1j*fy*L*Dy./(2*pi) )',0,length(ks)^2,length(ks)^2);
B=Jx+Jy;

BB=B([1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end],[1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end]);
KK=K([1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end],[1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end]);
%-------------------------------------------------------------------------
%         VV=-BB^-1(KK*P)   con VV(0,0)=0
%-------------------------------------------------------------------------
VV=zeros((length(ks)^2),1);
VV([1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end],1)...
    =-(BB\KK)*P(1,[1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end]).';
VV=reshape(VV,length(ks),length(ks));

[Vrec,~,~]=kspace2position1(VV,n,MM,0,L);
Vrec=Vrec-max(real(Vrec(:)));

MSE=immse(real(Vrec),real(Vorg));
Ns= [Ns (M-1)/2]; MSEs=[MSEs MSE];
P=reshape(Pkl,M^2,1);

Dy=repmat(n,1,length(n))';
Dx=sort(repmat(n,1,length(n)))';
k0k0=ceil(length(Jxnm)/2);
%-------------------------------------------------------------------------
%       Calculates (Fx,Fy) from guess value and the drift (Jx,Jy)
%-------------------------------------------------------------------------
Fx=real(Jxnm(k0k0,k0k0))+real(1i*2*pi*sum(Dx.*flipud(P).*VV(:)));
Fy=real(Jynm(k0k0,k0k0))+real(1i*2*pi*sum(Dy.*flipud(P).*VV(:))); % 2D
Fs=[Fx;Fy]; % 2D

diary('ResultsText.txt')
disp('------------------------')
disp('---- iterativo2D.m -----')
disp(['MSE=' num2str(MSE) ', K=' num2str((M-1)/2)])
disp(['fx=' num2str(fx,7) ', Fx=' num2str(Fx,7)])
disp(['fy=' num2str(fy,7) ', Fy=' num2str(Fy,7)]) 
disp('------------------------')
fileID = fopen('EvolutionIterations1.txt','at');
fprintf(fileID,'%2.12f %2.12f %2.12f %2.12f %2.12f\n',MSE,fx,Fx,fy,Fy);
fclose(fileID);
dlmwrite('Vrec.txt',real(Vrec));
end
