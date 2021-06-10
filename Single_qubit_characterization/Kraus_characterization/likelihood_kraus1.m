function [L] = likelihood_kraus1(M0z,M1z,M0x,M1x,M0y,M1y,POVM,rho,xm)
%Maximum log-likelihood functional for Kraus operators. 


%Obtain the Kraus operators
K1=[xm(1,1)+1i*xm(1,2) xm(1,3)+1i*xm(1,4);
    xm(1,5)+1i*xm(1,6) xm(1,7)+1i*xm(1,8)];

K2=[xm(1,9)+1i*xm(1,10) xm(1,11)+1i*xm(1,12);
    xm(1,13)+1i*xm(1,14) xm(1,15)+1i*xm(1,16)];

K3=[xm(1,17)+1i*xm(1,18) xm(1,19)+1i*xm(1,20);
    xm(1,21)+1i*xm(1,22) xm(1,23)+1i*xm(1,24)];

K4=[xm(1,25)+1i*xm(1,26) xm(1,27)+1i*xm(1,28);
    xm(1,29)+1i*xm(1,30) xm(1,31)+1i*xm(1,32)];


%The measurements in the z-basis 
f10=M0z(1,1);
f20=M0z(1,2);
f30=M0z(1,3);
f40=M0z(1,4);
f50=M0z(1,5);
f60=M0z(1,6);

f11=M1z(1,1);
f21=M1z(1,2);
f31=M1z(1,3);
f41=M1z(1,4);
f51=M1z(1,5);
f61=M1z(1,6);

c=1;

Lz=likelihood_kraus2(f10,f20,f30,f40,f50,f60,POVM,K1,K2,K3,K4,rho,c)+likelihood_kraus2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,K1,K2,K3,K4,rho,c);

%The measurements in the x-basis 
f10=M0x(1,1);
f20=M0x(1,2);
f30=M0x(1,3);
f40=M0x(1,4);
f50=M0x(1,5);
f60=M0x(1,6);

f11=M1x(1,1);
f21=M1x(1,2);
f31=M1x(1,3);
f41=M1x(1,4);
f51=M1x(1,5);
f61=M1x(1,6);

c=2;

Lx=likelihood_kraus2(f10,f20,f30,f40,f50,f60,POVM,K1,K2,K3,K4,rho,c)+likelihood_kraus2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,K1,K2,K3,K4,rho,c);

%The measurements in the y-basis 
f10=M0y(1,1);
f20=M0y(1,2);
f30=M0y(1,3);
f40=M0y(1,4);
f50=M0y(1,5);
f60=M0y(1,6);

f11=M1y(1,1);
f21=M1y(1,2);
f31=M1y(1,3);
f41=M1y(1,4);
f51=M1y(1,5);
f61=M1y(1,6);

c=3;

Ly=likelihood_kraus2(f10,f20,f30,f40,f50,f60,POVM,K1,K2,K3,K4,rho,c)+likelihood_kraus2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,K1,K2,K3,K4,rho,c);

%The -(log-likelihood)

L=-(Lz+Lx+Ly);

%Due to numerical imprecision Matlab can give imaginary output. 
%The fix below is to give warning about this if the imaginary error is too
%large and stear the optimization away from such points. Can be adjusted. 
if abs(imag(L))<10^(-8) 
L=real(L);
else
L=abs(L*100000);
disp('imaginary output');

end
    
end

