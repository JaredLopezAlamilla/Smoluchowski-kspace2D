function [KK,BB,P]=InversionMatrix(Pkl,M,fx,fy)
%-------------------------------------------------------------------------
%            BUILDS    SPARE DIAGONAL MATRIX FOR F-P EQ. 2D
%-------------------------------------------------------------------------
k=floor(M/2); ks=-k:1:k;
W=spdiags((ks.*(ks-0))',0,length(ks),length(ks));
D=spdiags(ones(1,length(ks))',0,length(ks),length(ks));

for l=1:k           
    Dg1=spdiags(ones(1,length(ks))',l,length(ks),length(ks));   
    Dg2=spdiags(ones(1,length(ks))',-l,length(ks),length(ks)); 
    Wg1=spdiags(circshift((ks.*(ks+l)),[0,l])',l,length(ks),length(ks));
    Wg2=spdiags(circshift((ks.*(ks-l)),[0,-l])',-l,length(ks),length(ks));
    
    D=D+Dg1+Dg2;      
    W=W+Wg1+Wg2;
end

D=repmat(D,length(ks));  
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

K=spdiags((Dx.*Dx-1j*fx*Dx./(2*pi)+Dy.*Dy-1j*fy*Dy./(2*pi) )',0,length(ks)^2,length(ks)^2);

B=Jx+Jy;

BB=B([1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end],[1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end]);
KK=K([1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end],[1:floor(length(ks)^2/2) ceil(length(ks)^2/2)+1:end]);
end
