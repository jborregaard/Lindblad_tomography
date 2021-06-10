%% Script to obtain the measurement SPAM errors
%In this script we use a maximum likelihood to extract the SPAM errors for two qubits 

%Load qubit data file giving measurement results for two qubits for t=0. (path is
%arbitrary)
load('Data/Data2020_2/Qb1/t_0.mat');
load('Data/Data2020_2/Qb2/t_0.mat');

U=1000; %Number of measurements per combination of input state and measurement basis. 

%In our case the measurement data was stored as arrays like m1(ii,jj,kk,ll,tt) where 
%ii=1,2,3,4,5,6 -> qubit 1 in |0>,|1>,|+>,|->,|+i>,|-i>
%jj=1,2,3,4,5,6 -> qubit 2 in |0>,|1>,|+>,|->,|+i>,|-i>
%kk=1,2,3 -> measurement of qubit 1 in basis x,y,z
%ll=1,2,3 -> measurement of qubit 2 in basis x,y,z
%tt=1:U -> the number of repetitions
%m1 is the measurement outcome of qubit 1 (0/1) and m2 is the measurement outcome
%of qubit 2 (0/1) 

%This measurement data is put into into m00, m01,m10,m11. Here mxy means outcome x for qubit 1
%and outcome y for qubit 2. The arrays are indexed with ii,jj,kk,ll with
%the same meaning as above. m00(1,1,1,1) gives the number of measurements with outcome 00 for input states |0>|0> and measurement
%in the z-basis for both qubits. 

m00=zeros(6,6,3,3);
m01=zeros(6,6,3,3);
m10=zeros(6,6,3,3);
m11=zeros(6,6,3,3);
for ii=1:6
    for jj=1:6
        for kk=1:3
            for ll=1:3
                for tt=1:U
                if m1(ii,jj,kk,ll,tt)==0 && m2(ii,jj,kk,ll,tt)==0
                    m00(ii,jj,kk,ll)=m00(ii,jj,kk,ll)+1;
                end
                if m1(ii,jj,kk,ll,tt)==0 && m2(ii,jj,kk,ll,tt)==1
                    m01(ii,jj,kk,ll)=m01(ii,jj,kk,ll)+1;
                end
                if m1(ii,jj,kk,ll,tt)==1 && m2(ii,jj,kk,ll,tt)==0
                    m10(ii,jj,kk,ll)=m10(ii,jj,kk,ll)+1;
                end
                if m1(ii,jj,kk,ll,tt)==1 && m2(ii,jj,kk,ll,tt)==1
                    m11(ii,jj,kk,ll)=m11(ii,jj,kk,ll)+1;
                end
                end
            end
        end
    end
end

%% Optimization

%Fmincon specifications such as constraints (see constraintpovm file)
nonlcon=@constraintpovm;
options = optimoptions('fmincon','MaxFunctionEvaluations',20000);

Sample=1000; %number of different starting points for optimization
XM=cell(1,Sample); %will store the output vector of optimization
L=zeros(1,Sample); %will store the output loglikelihood

for kk=1:Sample 
%Sample random POVM

%Put POVM elements into optimization vector for fmincon (using cholesky
%decomposition)
M000=ctranspose(chol(nearestSPD(M00)));
M001=ctranspose(chol(nearestSPD(M01)));
M010=ctranspose(chol(nearestSPD(M10)));
M011=ctranspose(chol(nearestSPD(M11)));

x00(1,1)=M000(1,1);
x00(1,2)=M000(2,2);
x00(1,3)=M000(3,3);
x00(1,4)=M000(4,4);
x00(1,5)=real(M000(2,1));
x00(1,6)=imag(M000(2,1));
x00(1,7)=real(M000(3,2));
x00(1,8)=imag(M000(3,2));
x00(1,9)=real(M000(4,3));
x00(1,10)=imag(M000(4,3));
x00(1,11)=real(M000(3,1));
x00(1,12)=imag(M000(3,1));
x00(1,13)=real(M000(4,2));
x00(1,14)=imag(M000(4,2));
x00(1,15)=real(M000(4,1));
x00(1,16)=imag(M000(4,1));


x01(1,1)=M001(1,1);
x01(1,2)=M001(2,2);
x01(1,3)=M001(3,3);
x01(1,4)=M001(4,4);
x01(1,5)=real(M001(2,1));
x01(1,6)=imag(M001(2,1));
x01(1,7)=real(M001(3,2));
x01(1,8)=imag(M001(3,2));
x01(1,9)=real(M001(4,3));
x01(1,10)=imag(M001(4,3));
x01(1,11)=real(M001(3,1));
x01(1,12)=imag(M001(3,1));
x01(1,13)=real(M001(4,2));
x01(1,14)=imag(M001(4,2));
x01(1,15)=real(M001(4,1));
x01(1,16)=imag(M001(4,1));

x10(1,1)=M010(1,1);
x10(1,2)=M010(2,2);
x10(1,3)=M010(3,3);
x10(1,4)=M010(4,4);
x10(1,5)=real(M010(2,1));
x10(1,6)=imag(M010(2,1));
x10(1,7)=real(M010(3,2));
x10(1,8)=imag(M010(3,2));
x10(1,9)=real(M010(4,3));
x10(1,10)=imag(M010(4,3));
x10(1,11)=real(M010(3,1));
x10(1,12)=imag(M010(3,1));
x10(1,13)=real(M010(4,2));
x10(1,14)=imag(M010(4,2));
x10(1,15)=real(M010(4,1));
x10(1,16)=imag(M010(4,1));

x11(1,1)=M011(1,1);
x11(1,2)=M011(2,2);
x11(1,3)=M011(3,3);
x11(1,4)=M011(4,4);
x11(1,5)=real(M011(2,1));
x11(1,6)=imag(M011(2,1));
x11(1,7)=real(M011(3,2));
x11(1,8)=imag(M011(3,2));
x11(1,9)=real(M011(4,3));
x11(1,10)=imag(M011(4,3));
x11(1,11)=real(M011(3,1));
x11(1,12)=imag(M011(3,1));
x11(1,13)=real(M011(4,2));
x11(1,14)=imag(M011(4,2));
x11(1,15)=real(M011(4,1));
x11(1,16)=imag(M011(4,1));


X0(1,1:16)=x00;
X0(1,17:32)=x01;
X0(1,33:48)=x10;
X0(1,49:64)=x11;

x0=real(X0);
    
%Sample random initial matrix

%Input density matrix to optimization vector (using cholesky decomposition)

rho0p=ctranspose(chol(rho0+1e-13*eye(size(rho0))));
x0(1,65)=rho0p(1,1);
x0(1,66)=rho0p(2,2);
x0(1,67)=rho0p(3,3);
x0(1,68)=rho0p(4,4);
x0(1,69)=real(rho0p(2,1));
x0(1,70)=imag(rho0p(2,1));
x0(1,71)=real(rho0p(3,2));
x0(1,72)=imag(rho0p(3,2));
x0(1,73)=real(rho0p(4,3));
x0(1,74)=imag(rho0p(4,3));
x0(1,75)=real(rho0p(3,1));
x0(1,76)=imag(rho0p(3,1));
x0(1,77)=real(rho0p(4,2));
x0(1,78)=imag(rho0p(4,2));
x0(1,79)=real(rho0p(4,1));
x0(1,80)=imag(rho0p(4,1)); 

x0(1,1)=rho0p(1,1);
x0(1,2)=rho0p(2,2);
x0(1,3)=rho0p(3,3);
x0(1,4)=rho0p(4,4);
x0(1,5)=real(rho0p(2,1));
x0(1,6)=imag(rho0p(2,1));
x0(1,7)=real(rho0p(3,2));
x0(1,8)=imag(rho0p(3,2));
x0(1,9)=real(rho0p(4,3));
x0(1,10)=imag(rho0p(4,3));
x0(1,11)=real(rho0p(3,1));
x0(1,12)=imag(rho0p(3,1));
x0(1,13)=real(rho0p(4,2));
x0(1,14)=imag(rho0p(4,2));
x0(1,15)=real(rho0p(4,1));
x0(1,16)=imag(rho0p(4,1));  
        
[xm,L0]=fmincon(@(x)likelihood_povm1(m00,m01,m10,m11,x),x0,[],[],[],[],[],[],nonlcon,options);
L(1,kk)=-L0;
XM{1,kk}=xm;
end

[Lmax,pos]=max(L); %find maximum loglikelihood.
nummaxL=find(L==Lmax);
if length(nummaxL)>1
disp('multiple minima'); %To give information about multiple maximas found. 
end
Xm=XM{1,pos}; %The optimal solution

%Extract the POVM elements from optim vector
Xm00=Xm(1,1:16);
Xm01=Xm(1,17:32);
Xm10=Xm(1,33:48);
Xm11=Xm(1,49:64);

%The POVM's
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

% Extracting the initial state from optim vector: 

Xmp=Xm(1,65:80);
Xmp=Xm;

%The density matrix

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

%Save output
save('SPAM2020_3.mat','M00','M01','M10','M11','rho','Lmax','-mat');