function [rho] = samplerho(epsilon)
%%Sample rho
%One can choose to look at deviation from a perfect state with various
%strength epsilon. 

Q = orth(randn(4) + randn(4)*1i);
D = diag(abs(randn(4, 1)));
A = Q*D*Q';
rhon=1/trace(A)*A;
rho0=nearestSPD(rhon);

rhop=zeros(4,4);
rhop(1,1)=1;
rho=(1-epsilon)*rhop+epsilon*rho0;
end

