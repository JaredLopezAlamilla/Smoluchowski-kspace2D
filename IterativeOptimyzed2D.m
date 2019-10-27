Ns=[];MSEs=[];
epsilon=[0.0007;0.0007]; % 2D case
options = optimset('TolFun',epsilon);

fileID = fopen('EvolutionIterations1.txt','at');
fprintf(fileID, '%s\t\t%s\t\t%s\t\t%s\t\t%s\n','MSEpot','fx','Fx','fy','Fy');
fclose(fileID);

[Fs,fval,exitflag,output]= fsolve(@(F) abs(F-iterativo2D(F,P_pos,V2D,Jxkq,Jykq,Ns,MSEs)),[0;0],options); % 2D case
