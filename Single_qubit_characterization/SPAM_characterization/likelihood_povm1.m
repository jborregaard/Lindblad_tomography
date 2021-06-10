function [L] = likelihood_povm1(M0z,M1z,M0x,M1x,M0y,M1y,xm)
%Maximum likelihood functional for extraxting SPAM errors.

%The POVM elements

POVMo(1,1)=xm(1,1);
POVMo(2,1)=xm(1,2)+1i*xm(1,3);
POVMo(2,2)=xm(1,4);

POVM=POVMo*ctranspose(POVMo); 

%The initial state
rhoo(1,1)=xm(1,5);
rhoo(2,1)=xm(1,7)+1i*xm(1,8);
rhoo(2,2)=xm(1,6);

rhop=rhoo*ctranspose(rhoo);
rho=1/trace(rhop)*rhop;

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

Lz=likelihood_povm2(f10,f20,f30,f40,f50,f60,POVM,rho,c)+likelihood_povm2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,rho,c);

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

Lx=likelihood_povm2(f10,f20,f30,f40,f50,f60,POVM,rho,c)+likelihood_povm2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,rho,c);

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

Ly=likelihood_povm2(f10,f20,f30,f40,f50,f60,POVM,rho,c)+likelihood_povm2(f11,f21,f31,f41,f51,f61,eye(2)-POVM,rho,c);

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

