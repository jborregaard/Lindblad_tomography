function [L] = likelihood_povm1(m00,m01,m10,m11,X)
%Maximum likelihood functional for the SPAM errors for 2 qubits

Xm=X;

Xm00=Xm(1,1:16);
Xm01=Xm(1,17:32);
Xm10=Xm(1,33:48);
Xm11=Xm(1,49:64);

Xmp=Xm(1,65:80);


%The POVM elements

L00(1,1)=abs(Xm00(1,1));
L00(2,2)=abs(Xm00(1,2));
L00(3,3)=abs(Xm00(1,3));
L00(4,4)=abs(Xm00(1,4));
L00(2,1)=Xm00(1,5)+1i*Xm00(1,6);
L00(3,2)=Xm00(1,7)+1i*Xm00(1,8);
L00(4,3)=Xm00(1,9)+1i*Xm00(1,10);
L00(3,1)=Xm00(1,11)+1i*Xm00(1,12);
L00(4,2)=Xm00(1,13)+1i*Xm00(1,14);
L00(4,1)=Xm00(1,15)+1i*Xm00(1,16);

M00=L00*ctranspose(L00);

L01(1,1)=abs(Xm01(1,1));
L01(2,2)=abs(Xm01(1,2));
L01(3,3)=abs(Xm01(1,3));
L01(4,4)=abs(Xm01(1,4));
L01(2,1)=Xm01(1,5)+1i*Xm01(1,6);
L01(3,2)=Xm01(1,7)+1i*Xm01(1,8);
L01(4,3)=Xm01(1,9)+1i*Xm01(1,10);
L01(3,1)=Xm01(1,11)+1i*Xm01(1,12);
L01(4,2)=Xm01(1,13)+1i*Xm01(1,14);
L01(4,1)=Xm01(1,15)+1i*Xm01(1,16);

M01=L01*ctranspose(L01);

L10(1,1)=abs(Xm10(1,1));
L10(2,2)=abs(Xm10(1,2));
L10(3,3)=abs(Xm10(1,3));
L10(4,4)=abs(Xm10(1,4));
L10(2,1)=Xm10(1,5)+1i*Xm10(1,6);
L10(3,2)=Xm10(1,7)+1i*Xm10(1,8);
L10(4,3)=Xm10(1,9)+1i*Xm10(1,10);
L10(3,1)=Xm10(1,11)+1i*Xm10(1,12);
L10(4,2)=Xm10(1,13)+1i*Xm10(1,14);
L10(4,1)=Xm10(1,15)+1i*Xm10(1,16);

M10=L10*ctranspose(L10);

L11(1,1)=abs(Xm11(1,1));
L11(2,2)=abs(Xm11(1,2));
L11(3,3)=abs(Xm11(1,3));
L11(4,4)=abs(Xm11(1,4));
L11(2,1)=Xm11(1,5)+1i*Xm11(1,6);
L11(3,2)=Xm11(1,7)+1i*Xm11(1,8);
L11(4,3)=Xm11(1,9)+1i*Xm11(1,10);
L11(3,1)=Xm11(1,11)+1i*Xm11(1,12);
L11(4,2)=Xm11(1,13)+1i*Xm11(1,14);
L11(4,1)=Xm11(1,15)+1i*Xm11(1,16);

M11=L11*ctranspose(L11);

%The initial density matrix

rhop(1,1)=abs(Xmp(1,1));
rhop(2,2)=abs(Xmp(1,2));
rhop(3,3)=abs(Xmp(1,3));
rhop(4,4)=abs(Xmp(1,4));
rhop(2,1)=Xmp(1,5)+1i*Xmp(1,6);
rhop(3,2)=Xmp(1,7)+1i*Xmp(1,8);
rhop(4,3)=Xmp(1,9)+1i*Xmp(1,10);
rhop(3,1)=Xmp(1,11)+1i*Xmp(1,12);
rhop(4,2)=Xmp(1,13)+1i*Xmp(1,14);
rhop(4,1)=Xmp(1,15)+1i*Xmp(1,16);

rho1=rhop*ctranspose(rhop);

rho=1/trace(rho1)*rho1;

%The measurements in the zz-basis 
ba1=1; %specifies measurement basis of qubit 1 (same below)
ba2=1; %specifies measurement basis of qubit 2 (same below)
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=1; %specifies measurement basis for likelihood_povm2 (same below)

Lzz00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzz01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzz10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzz11=likelihood_povm2(f,M11,rho,c);

Lzz=Lzz00+Lzz01+Lzz10+Lzz11;

%The measurements in the zx-basis 
ba1=1;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=2;

Lzx00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzx01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzx10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzx11=likelihood_povm2(f,M11,rho,c);

Lzx=Lzx00+Lzx01+Lzx10+Lzx11;

%The measurements in the zy-basis 
ba1=1;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=3;

Lzy00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzy01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzy10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzy11=likelihood_povm2(f,M11,rho,c);

Lzy=Lzy00+Lzy01+Lzy10+Lzy11;
%The measurements in the xz-basis 
ba1=3;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=4;

Lxz00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxz01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxz10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxz11=likelihood_povm2(f,M11,rho,c);

Lxz=Lxz00+Lxz01+Lxz10+Lxz11;

%The measurements in the xx-basis 
ba1=3;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=5;

Lxx00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxx01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxx10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxx11=likelihood_povm2(f,M11,rho,c);

Lxx=Lxx00+Lxx01+Lxx10+Lxx11;

%The measurements in the xy-basis 
ba1=3;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=6;

Lxy00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxy01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxy10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxy11=likelihood_povm2(f,M11,rho,c);

Lxy=Lxy00+Lxy01+Lxy10+Lxy11;

%The measurements in the yz-basis 
ba1=2;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=7;

Lyz00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyz01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyz10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyz11=likelihood_povm2(f,M11,rho,c);

Lyz=Lyz00+Lyz01+Lyz10+Lyz11;

%The measurements in the yx-basis 
ba1=2;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=8;

Lyx00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyx01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyx10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyx11=likelihood_povm2(f,M11,rho,c);

Lyx=Lyx00+Lyx01+Lyx10+Lyx11;

%The measurements in the yy-basis 
ba1=2;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=9;

Lyy00=likelihood_povm2(f,M00,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyy01=likelihood_povm2(f,M01,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyy10=likelihood_povm2(f,M10,rho,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyy11=likelihood_povm2(f,M11,rho,c);

Lyy=Lyy00+Lyy01+Lyy10+Lyy11;

%The -(log-likelihood) 

L=-(Lzz+Lzy+Lzx+Lxz+Lxy+Lxx+Lyz+Lyy+Lyx); 

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

