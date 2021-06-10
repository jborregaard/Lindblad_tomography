function [Out] = likelihood_Lindblad1(m00t,m01t,m10t,m11t,M00,M01,M10,M11,S,rho,X)
%Maximum likelihood functional for the Lindbladian. 

%U is the number of time steps in the measurement data
U=160;

%Getting the Lindblad matrix from the optim vector
X2=zeros(1,241);
X2(1,1:225)=X(1,17:241);
X2(1,226:241)=X(1,1:16);

xm=X2;
odd=[1:2:31];
hnmp=zeros(15,15);
for ttt=1:15
for jjj=ttt:15
if ttt==jjj
hnmp(jjj,ttt)=abs(xm(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+1));
else
hnmp(jjj,ttt)=xm(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+2*(jjj-ttt))+1i*xm(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+1+2*(jjj-ttt));
end
end
end

hnm=hnmp*ctranspose(hnmp);

%Geting the Hamiltonian from the optim vector
H(1,1)=xm(1,226);
H(1,2)=xm(1,227)+1i*xm(1,228);
H(2,1)=conj(H(1,2));
H(1,3)=xm(1,229)+1i*xm(1,230);
H(3,1)=conj(H(1,3));
H(1,4)=xm(1,231)+1i*xm(1,232);
H(4,1)=conj(H(1,4));
H(2,2)=xm(1,233);
H(2,3)=xm(1,234)+1i*xm(1,235);
H(3,2)=conj(H(3,2));
H(2,4)=xm(1,236)+1i*xm(1,237);
H(4,2)=conj(H(2,4));
H(3,3)=xm(1,238);
H(3,4)=xm(1,239)+1i*xm(1,240);
H(4,3)=conj(H(3,4));
H(4,4)=xm(1,241);


%the single qubit Gates for initial state preparation:
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

y=cell(U,36);
y{1,1}=p00;
y{1,2}=p01;
y{1,3}=p0p;
y{1,4}=p0m;
y{1,5}=p0ip;
y{1,6}=p0im;
y{1,7}=p10;
y{1,8}=p11;
y{1,9}=p1p;
y{1,10}=p1m;
y{1,11}=p1ip;
y{1,12}=p1im;
y{1,13}=pp0;
y{1,14}=pp1;
y{1,15}=ppp;
y{1,16}=ppm;
y{1,17}=ppip;
y{1,18}=ppim;
y{1,19}=pm0;
y{1,20}=pm1;
y{1,21}=pmp;
y{1,22}=pmm;
y{1,23}=pmip;
y{1,24}=pmim;
y{1,25}=pip0;
y{1,26}=pip1;
y{1,27}=pipp;
y{1,28}=pipm;
y{1,29}=pipip;
y{1,30}=pipim;
y{1,31}=pim0;
y{1,32}=pim1;
y{1,33}=pimp;
y{1,34}=pimm;
y{1,35}=pimip;
y{1,36}=pimim;

%Load the data time vector and use unit of 1 microsec.
load('../Data/Data2020_2/time.mat');
time=time*10^(6);
Loglike=zeros(1,U);

%Construct liouvillian superoperator from the Lindblad matrix
Id=eye(4);
L=zeros(16,16);
for ii=1:15
for jj=1:15
Lp=hnm(ii,jj).*(kron(S{1,ii},Id)*kron(Id,transpose(ctranspose(S{1,jj})))-1/2*(kron(ctranspose(S{1,jj})*S{1,ii},Id)+kron(Id,transpose(ctranspose(S{1,jj})*S{1,ii}))));    
L=L+Lp;
end
end

% Add the Hamiltonian evolution
L=1i*(kron(Id,transpose(H))-kron(H,Id))+L;

%Get the density matrices after evolution under the Louvillian
rhop=cell(1,36);
for kk=1:36
Yp=zeros(length(time),16);
for nn=1:length(time)
Yp(nn,:)=expm(L*time(1,nn))*reshape(y{1,kk},16,1); 
end
rhop{1,kk}=Yp;
end

%Calculate the log-likelihood for each timestep. 
for ii=1:U       
for kk=1:36
Ypp=rhop{1,kk};    
y{ii+1,kk}=reshape(Ypp(ii+1,:),4,4);
end

p00f=y{ii+1,1};
p01f=y{ii+1,2};
p0pf=y{ii+1,3};
p0mf=y{ii+1,4};
p0ipf=y{ii+1,5};
p0imf=y{ii+1,6};

p10f=y{ii+1,7};
p11f=y{ii+1,8};
p1pf=y{ii+1,9};
p1mf=y{ii+1,10};
p1ipf=y{ii+1,11};
p1imf=y{ii+1,12};

pp0f=y{ii+1,13};
pp1f=y{ii+1,14};
pppf=y{ii+1,15};
ppmf=y{ii+1,16};
ppipf=y{ii+1,17};
ppimf=y{ii+1,18};

pm0f=y{ii+1,19};
pm1f=y{ii+1,20};
pmpf=y{ii+1,21};
pmmf=y{ii+1,22};
pmipf=y{ii+1,23};
pmimf=y{ii+1,24};

pip0f=y{ii+1,25};
pip1f=y{ii+1,26};
pippf=y{ii+1,27};
pipmf=y{ii+1,28};
pipipf=y{ii+1,29};
pipimf=y{ii+1,30};

pim0f=y{ii+1,31};
pim1f=y{ii+1,32};
pimpf=y{ii+1,33};
pimmf=y{ii+1,34};
pimipf=y{ii+1,35};
pimimf=y{ii+1,36};

m00=m00t{1,ii};
m01=m01t{1,ii};
m10=m10t{1,ii};
m11=m11t{1,ii};

%The measurements in the zz-basis 
ba1=1;%specifies measurement basis of qubit 1 (same below)
ba2=1;%specifies measurement basis of qubit 2 (same below)
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=1; %specifies measurement basis for likelihood_povm2 (same below)

Lzz00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzz01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzz10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzz11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzz=Lzz00+Lzz01+Lzz10+Lzz11;

%The measurements in the zx-basis 
ba1=1;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=2;

Lzx00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzx01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzx10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzx11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzx=Lzx00+Lzx01+Lzx10+Lzx11;

%The measurements in the zy-basis 
ba1=1;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=3;

Lzy00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lzy01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lzy10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lzy11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lzy=Lzy00+Lzy01+Lzy10+Lzy11;
%The measurements in the xz-basis 
ba1=3;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=4;

Lxz00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxz01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxz10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxz11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxz=Lxz00+Lxz01+Lxz10+Lxz11;

%The measurements in the xx-basis 
ba1=3;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=5;

Lxx00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxx01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxx10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxx11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxx=Lxx00+Lxx01+Lxx10+Lxx11;

%The measurements in the xy-basis 
ba1=3;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=6;

Lxy00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lxy01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lxy10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lxy11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lxy=Lxy00+Lxy01+Lxy10+Lxy11;

%The measurements in the yz-basis 
ba1=2;
ba2=1;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=7;

Lyz00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyz01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyz10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyz11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lyz=Lyz00+Lyz01+Lyz10+Lyz11;

%The measurements in the yx-basis 
ba1=2;
ba2=3;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=8;

Lyx00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyx01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyx10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyx11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lyx=Lyx00+Lyx01+Lyx10+Lyx11;

%The measurements in the yy-basis 
ba1=2;
ba2=2;
f=zeros(6,6);
f(1:6,1:6)=reshape(m00(1:6,1:6,ba1,ba2),6,6);
c=9;

Lyy00=likelihood_Lindblad2(f,M00,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m01(1:6,1:6,ba1,ba2),6,6);

Lyy01=likelihood_Lindblad2(f,M01,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m10(1:6,1:6,ba1,ba2),6,6);

Lyy10=likelihood_Lindblad2(f,M10,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

f=zeros(6,6);
f(1:6,1:6)=reshape(m11(1:6,1:6,ba1,ba2),6,6);

Lyy11=likelihood_Lindblad2(f,M11,p00f,p01f,p0pf,p0mf,p0ipf,p0imf,p10f,p11f,p1pf,p1mf,p1ipf,p1imf,pp0f,pp1f,pppf,ppmf,ppipf,ppimf,pm0f,pm1f,pmpf,pmmf,pmipf,pmimf,pip0f,pip1f,pippf,pipmf,pipipf,pipimf,pim0f,pim1f,pimpf,pimmf,pimipf,pimimf,c);

Lyy=Lyy00+Lyy01+Lyy10+Lyy11;

%The -(log-likelihood)
Loglike(1,ii)=-(Lzz+Lzy+Lzx+Lxz+Lxy+Lxx+Lyz+Lyy+Lyx); 

%Due to numerical imprecision Matlab can give imaginary output. 
%The fix below is to give warning about this if the imaginary error is too
%large and stear the optimization away from such points. Can be adjusted. 
if abs(imag(Loglike(1,ii)))/abs(real(Loglike(1,ii)))<10^(-6)     
Loglike(1,ii)=real(Loglike(1,ii));
else
disp('Imaginary part - numerical precision problem');  
Loglike(1,ii)=10^10*abs(Loglike(1,ii));
end
    
end
%The total -(log-likelihood)
Out=sum(Loglike);   
end

