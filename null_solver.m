function [Pnm,Jxnm,Jynm,vx,vy,vorticity,Norm]=null_solver(A,Jxk,Jyk,N,L,n)
% Solves the nullspace of A and returns P_nm matrix

    [Q,~,R]=qr(sparse(A.'),0);  % Creates upper triangular matrix such A = Q*R
    
    null_vec=conj(Q(:,end)); % Conjugated in vector shape
   
    midNvec=round((N^2)/2);
    Norm=null_vec(midNvec,1)/L^2;
    null_vec=null_vec/null_vec(midNvec,1)/L^2; % Normalise
    
    Pnm=reshape(null_vec,[N N]); % This creates Pnm
    
    Jxnm=reshape(-1i*2*pi*(Jxk*null_vec),[N N]); 
    Jynm=reshape(-1i*2*pi*(Jyk*null_vec),[N N]);
    vorticity=1i*(-repmat(n,N,1).*Jxnm + repmat(n,N,1)'.*Jynm);
    mid=round(length(Jxnm)/2);
    vx=real(Jxnm(mid,mid))*L; 
    vy=real(Jynm(mid,mid))*L;
end
