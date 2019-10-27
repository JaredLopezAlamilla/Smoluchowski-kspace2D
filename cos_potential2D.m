function [V0] = cos_potential2D(x,y,Vx,Vy,Vz)
% Returns the values of the cosine potential
% New custom periodic potentials can be created in a similar way 
    V0 =Vx*cos(2*pi*x)+Vy*cos(2*pi*y)+Vz*cos(2*pi*(x-0.5-y)); 
end
