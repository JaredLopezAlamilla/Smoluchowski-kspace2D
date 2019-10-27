function [pos_matrix,X,Y] = kspace2position1(k_matrix,n,res,periods,T)
% Convert matrix from k-space to position-space    

    x=T*(0:res-1)/res;
    
    for k=1:periods-1
        x=[x k+T*(0:res-1)/res];
    end
    
    [X,Y]= meshgrid(x);
    expnx = exp(1i*(2*pi/T)*n'*x);
    expny = exp(1i*(2*pi/T)*n'*x);
    pos_matrix = ((k_matrix*expnx).'*expny).';

end
