function [POVM] = samplepovm(epsilon)
%%Sample rho
%One can choose to look at deviation from a perfect POVM with various
%strength epsilon. 
loopCounter=1;
g=1;
while g==1
Q = orth(randn(4) + randn(4)*1i);
D = diag(abs(randn(4, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

rhoi=zeros(4,4);
rhoi(1,1)=1;
M00=(1-epsilon)*rhoi+epsilon*rho0;

Q = orth(randn(4) + randn(4)*1i);
D = diag(abs(randn(4, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

rhoi=zeros(4,4);
rhoi(2,2)=1;
M01=(1-epsilon)*rhoi+epsilon*rho0;

Q = orth(randn(4) + randn(4)*1i);
D = diag(abs(randn(4, 1)));
A = Q*D*Q';
rho0n=nearestSPD(A);
rho0=1/trace(rho0n)*rho0n;

rhoi=zeros(4,4);
rhoi(3,3)=1;
M10=(1-epsilon)*rhoi+epsilon*rho0;

E=eig(eye(4)-M00-M01-M10);
if E(1,1)>=0 && E(2,1)>=0 && E(3,1)>=0 && E(4,1)>=0
g=0;
M11=eye(4)-M00-M01-M10;
end

if loopCounter >= 3000
disp('too many rounds');    
    break;
end
  loopCounter  = loopCounter  + 1;  
end
POVM{1,1}=M00;
POVM{1,2}=M01;
POVM{1,3}=M10;
POVM{1,4}=M11;

end

