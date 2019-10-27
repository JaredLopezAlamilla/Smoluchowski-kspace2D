function [A,Jxk,Jyk]=find_matrix(Vkl,n,M,kBTx,kBTy,gx,gy,fx,fy,L)
%-------------------------------------------------------------------------
%         Creates the matrix  A*Pnm=0 with k-space potential Vkl.
%-------------------------------------------------------------------------
    m=n;
    
   [Jx0,Jy0]=deal(zeros(M^2));  % creates a zeros MxM matrix for Jx0 and Jy0          
    Vkl_vec=Vkl(:); % compres to a single vector
    
    k_vec=sort(repmat(n,1,M),2)'; % generate a vector of size Mx1 with all entires n     
    l_vec=repmat(sort(m,2),1,M)'; % generate a vector of size Mx1 with all entires m
    
    R=-floor(M^2/2):floor(M^2/2); % generate a vector of size Mx1 with all entires -m/2 ... m/2
    
    for q=1:M^2 % re-orders Vkl so can compute Jx0 and Jy0 correctly        
      Vkl_shifted=circshift(Vkl_vec,[R(q) 0]);  
      k_shifted=circshift(k_vec,[R(q) 0]);
      l_shifted=circshift(l_vec,[R(q) 0]);
      Jx0(:,q)=Vkl_shifted.*k_shifted;
      Jy0(:,q)=Vkl_shifted.*l_shifted;
    end
    
    diagonal_x=kron(diag(kBTx*n + 1i*fx*L/2/pi),spdiags(ones(M,1),0,M,M)); 
    diagonal_y=kron(spdiags(ones(M,1),0,M,M),diag(kBTy*m + 1i*fy*L/2/pi));
    
    Jx1=(Jx0+diagonal_x)/gx;    
    Jy1=(Jy0+diagonal_y)/gy;
    
    mat_fix=zeros(M,M);
    
    for w=n
      mat_fix = mat_fix + diag(ones(1,M-abs(w)),w);
    end
    
    mat_fix2=kron(mat_fix,mat_fix);  
    
    Jxk=Jx1.*mat_fix2/L;
    Jyk=Jy1.*mat_fix2/L;
%-------------------------------------------------------------------------
%         A = i*n*Jxk + i*m*Jyk  matrix for FP eq in k-space
%-------------------------------------------------------------------------   
    Ax=sort((repmat(n,M^2,M)')).*Jxk; 
    Ay=(repmat(m,M^2,M))'.*Jyk; 
    A=Ax+Ay;  % 
end
