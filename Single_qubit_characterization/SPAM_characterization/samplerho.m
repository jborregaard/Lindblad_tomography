function [rho] = samplerho(epsilon)
%%Sample rho
%One can choose to look at deviation from a perfect state with various
%strength epsilon. 

Q = orth(randn(2) + randn(2)*1i);
D = diag(abs(randn(2, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

rho=(1-epsilon)*[1 0; 0 0]+epsilon*rho0;
end

