function [rhof] = evolution(K1,K2,K3,K4,rho)
%Function to get the evolution of rho given by the Kraus operators. 

rhof=K1*rho*ctranspose(K1)+K2*rho*ctranspose(K2)+K3*rho*ctranspose(K3)+K4*rho*ctranspose(K4);

end

