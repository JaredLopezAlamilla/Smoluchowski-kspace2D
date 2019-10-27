# Smoluchowski-kSpace Iterative

This repository contains a collection of matlab scripts for  solving  various physical problems
invloving solutions to the Smoluchowski equation. 

The smoluchowski equation is represented in k-Space and standard linear algebra methods
are applied to solve it.
## Calculate steady-state probability  using the energy landscape
        
Functions required: SteadyState2D_V2P.m, find_matrix.m, null_slover. m, kspace2position1.m  
### Use example:
`>> SteadyState2D_V2P;`

plotting the results

`>> figure(); contour(Xfine,Yfine,V2D);title('$V_0(x,y)$');`

`>> figure(); contour(Xfine,Yfine,P_pos);title('$P(x,y)$');`

## Reconstruct the energy landscape  using the steady-state distribution  and drift velocities 
Functions required:  IterativeOptimyzed2D.m, iterativo2D.m, InversionMatrix.m, SolverP2V.m, kspace2position1.m
### Use example:
`>> IterativeOptimyzed2D;`

`>>Vrec=dlmread('Vrec.txt');`

plotting the results
        
`>> figure(); contour(Xfine,Yfine,Vrec);title('$V_0^{\rm rec}(x,y)$');`
### comparing the model vs reconstruction
`>> MSEpot=immse(Vrec,V2D);`

   *Additionally you need the functions to generate the different model potentials or build them inside the file SteadyState2D_V2P.m  
