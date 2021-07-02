function [POVM] = samplepovm(epsilon)
%%Sample rho
%One can choose to look at deviation from a perfect POVM with various
%strength epsilon. 
loopCounter=1;
g=1;
while g==1
Q = orth(randn(2) + randn(2)*1i);
D = diag(abs(randn(2, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

POVM=(1-epsilon)*[1 0; 0 0]+epsilon*rho0;
E=eig(eye(2)-POVM);
if E(1,1)>=0 && E(2,1)>=0
g=0;
end

if loopCounter >= 3000
disp('too many rounds');    
    break;
end
  loopCounter  = loopCounter  + 1;  
end

end

