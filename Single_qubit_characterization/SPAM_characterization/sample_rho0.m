function [Out] = sample_rho0(~)
%Function that samples a random singel qubit density matrix 
%%Sample rho
Q = orth(randn(2) + randn(2)*1i);
D = diag(abs(randn(2, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

Out=rho0;
end

