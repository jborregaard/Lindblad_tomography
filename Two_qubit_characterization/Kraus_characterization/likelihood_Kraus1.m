function [L] = likelihood_Kraus1(m00,m01,m10,m11,M00,M01,M10,M11,rho,X)
%Maximum likelihood functional for the Kraus operators.



%the single qubit Gates for qubit preparation:
sigmax=[0 -1i;-1i 0];
hadp=1/sqrt(2)*[1 -1;1 1];
hadip=1/sqrt(2)*[1 1i; 1i 1];
hadm=1/sqrt(2)*[1 -1;1 1];
hadim=1/sqrt(2)*[1 -1i; -1i 1];


%The initial density matrices
%p00=|0>|0>
%p0p=|0>|+>
%p0m=|0>|->
%p0ip=|0>|+i>
%p0im=|0>|-i> etc. 

p00=rho;
p01=kron(eye(2),sigmax)*rho*ctranspose(kron(eye(2),sigmax));
p0p=kron(eye(2),hadp)*rho*ctranspose(kron(eye(2),hadp));
p0m=kron(eye(2),hadm)*rho*ctranspose(kron(eye(2),hadm));
p0ip=kron(eye(2),hadip)*rho*ctranspose(kron(eye(2),hadip));
p0im=kron(eye(2),hadim)*rho*ctranspose(kron(eye(2),hadim));

p10=kron(sigmax,eye(2))*rho*ctranspose(kron(sigmax,eye(2)));
p11=kron(sigmax,sigmax)*rho*ctranspose(kron(sigmax,sigmax));
p1p=kron(sigmax,hadp)*rho*ctranspose(kron(sigmax,hadp));
p1m=kron(sigmax,hadm)*rho*ctranspose(kron(sigmax,hadm));
p1ip=kron(sigmax,hadip)*rho*ctranspose(kron(sigmax,hadip));
p1im=kron(sigmax,hadim)*rho*ctranspose(kron(sigmax,hadim));

pp0=kron(hadp,eye(2))*rho*ctranspose(kron(hadp,eye(2)));
pp1=kron(hadp,sigmax)*rho*ctranspose(kron(hadp,sigmax));
ppp=kron(hadp,hadp)*rho*ctranspose(kron(hadp,hadp));
ppm=kron(hadp,hadm)*rho*ctranspose(kron(hadp,hadm));
ppip=kron(hadp,hadip)*rho*ctranspose(kron(hadp,hadip));
ppim=kron(hadp,hadim)*rho*ctranspose(kron(hadp,hadim));

pm0=kron(hadm,eye(2))*rho*ctranspose(kron(hadm,eye(2)));
pm1=kron(hadm,sigmax)*rho*ctranspose(kron(hadm,sigmax));
pmp=kron(hadm,hadp)*rho*ctranspose(kron(hadm,hadp));
pmm=kron(hadm,hadm)*rho*ctranspose(kron(hadm,hadm));
pmip=kron(hadm,hadip)*rho*ctranspose(kron(hadm,hadip));
pmim=kron(hadm,hadim)*rho*ctranspose(kron(hadm,hadim));

pip0=kron(hadip,eye(2))*rho*ctranspose(kron(hadip,eye(2)));
pip1=kron(hadip,sigmax)*rho*ctranspose(kron(hadip,sigmax));
pipp=kron(hadip,hadp)*rho*ctranspose(kron(hadip,hadp));
pipm=kron(hadip,hadm)*rho*ctranspose(kron(hadip,hadm));
pipip=kron(hadip,hadip)*rho*ctranspose(kron(hadip,hadip));
pipim=kron(hadip,hadim)*rho*ctranspose(kron(hadip,hadim));

pim0=kron(hadim,eye(2))*rho*ctranspose(kron(hadim,eye(2)));
pim1=kron(hadim,sigmax)*rho*ctranspose(kron(hadim,sigmax));
pimp=kron(hadim,hadp)*rho*ctranspose(kron(hadim,hadp));
pimm=kron(hadim,hadm)*rho*ctranspose(kron(hadim,hadm));
pimip=kron(hadim,hadip)*rho*ctranspose(kron(hadim,hadip));
pimim=kron(hadim,hadim)*rho*ctranspose(kron(hadim,hadim));

%The Kraus operators from the optim vector
K1=reshape(X(1,1:16),4,4)+1i*reshape(X(1,256+1:256+16),4,4);
K2=reshape(X(1,1+16:2*16),4,4)+1i*reshape(X(1,256+1+16:256+2*16),4,4);
K3=reshape(X(1,1+2*16:3*16),4,4)+1i*reshape(X(1,256+1+2*16:256+3*16),4,4);
K4=reshape(X(1,1+3*16:4*16),4,4)+1i*reshape(X(1,256+1+3*16:256+4*16),4,4);
K5=reshape(X(1,1+4*16:5*16),4,4)+1i*reshape(X(1,256+1+4*16:256+5*16),4,4);
K6=reshape(X(1,1+5*16:6*16),4,4)+1i*reshape(X(1,256+1+5*16:256+6*16),4,4);
K7=reshape(X(1,1+6*16:7*16),4,4)+1i*reshape(X(1,256+1+6*16:256+7*16),4,4);
K8=reshape(X(1,1+7*16:8*16),4,4)+1i*reshape(X(1,256+1+7*16:256+8*16),4,4);
K9=reshape(X(1,1+8*16:9*16),4,4)+1i*reshape(X(1,256+1+8*16:256+9*16),4,4);
K10=reshape(X(1,1+9*16:10*16),4,4)+1i*reshape(X(1,256+1+9*16:256+10*16),4,4);
K11=reshape(X(1,1+10*16:11*16),4,4)+1i*reshape(X(1,256+1+10*16:256+11*16),4,4);
K12=reshape(X(1,1+11*16:12*16),4,4)+1i*reshape(X(1,256+1+11*16:256+12*16),4,4);
K13=reshape(X(1,1+12*16:13*16),4,4)+1i*reshape(X(1,256+1+12*16:256+13*16),4,4);
K14=reshape(X(1,1+13*16:14*16),4,4)+1i*reshape(X(1,256+1+13*16:256+14*16),4,4);
K15=reshape(X(1,1+14*16:15*16),4,4)+1i*reshape(X(1,256+1+14*16:256+15*16),4,4);
K16=reshape(X(1,1+15*16:16*16),4,4)+1i*reshape(X(1,256+1+15*16:256+16*16),4,4);


%Get the evolution of the different starting states according to the Kraus
%operators
p00f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p00);
p01f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p01);
p0pf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p0p);
p0mf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p0m);
p0ipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p0ip);
p0imf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p0im);

p10f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p10);
p11f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p11);
p1pf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p1p);
p1mf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p1m);
p1ipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p1ip);
p1imf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,p1im);

pp0f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pp0);
pp1f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pp1);
pppf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,ppp);
ppmf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,ppm);
ppipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,ppip);
ppimf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,ppim);

pm0f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pm0);
pm1f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pm1);
pmpf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pmp);
pmmf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pmm);
pmipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pmip);
pmimf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pmim);

pip0f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pip0);
pip1f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pip1);
pippf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pipp);
pipmf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pipm);
pipipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pipip);
pipimf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pipim);

pim0f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pim0);
pim1f=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pim1);
pimpf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pimp);
pimmf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pimm);
pimipf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pimip);
pimimf=evolution(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,pimim);

%% The log likelihood functions
%The measurements in the zz-basis 
ba1=1; %specifies measurement basis of qubit 1 (same below)
ba2=1; %specifies measurement basis of qubit 2 (same below)
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=1; %specifies measurement basis for likelihood_povm2 (same below)

Lzz00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzz01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzz10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzz11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzz=Lzz00+Lzz01+Lzz10+Lzz11;

%The measurements in the zx-basis 
ba1=1;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=2;

Lzx00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzx01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzx10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzx11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzx=Lzx00+Lzx01+Lzx10+Lzx11;

%The measurements in the zy-basis 
ba1=1;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=3;

Lzy00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzy01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzy10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzy11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzy=Lzy00+Lzy01+Lzy10+Lzy11;
%The measurements in the xz-basis 
ba1=3;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=4;

Lxz00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxz01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxz10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxz11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxz=Lxz00+Lxz01+Lxz10+Lxz11;

%The measurements in the xx-basis 
ba1=3;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=5;

Lxx00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxx01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxx10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxx11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxx=Lxx00+Lxx01+Lxx10+Lxx11;

%The measurements in the xy-basis 
ba1=3;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=6;

Lxy00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxy01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxy10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxy11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxy=Lxy00+Lxy01+Lxy10+Lxy11;

%The measurements in the yz-basis 
ba1=2;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=7;

Lyz00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyz01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyz10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyz11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lyz=Lyz00+Lyz01+Lyz10+Lyz11;

%The measurements in the yx-basis 
ba1=2;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=8;

Lyx00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyx01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyx10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyx11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lyx=Lyx00+Lyx01+Lyx10+Lyx11;

%The measurements in the yy-basis 
ba1=2;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=9;

Lyy00=likelihood_Kraus2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyy01=likelihood_Kraus2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyy10=likelihood_Kraus2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyy11=likelihood_Kraus2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

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

