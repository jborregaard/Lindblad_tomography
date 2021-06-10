%% M-file to fit the Lindblad operators
%In this script we use a maximum likelihood to extract the Lindblad operators 
tic;
%Load extracted SPAM
load('../SPAM2020_1_0.mat');

%Ut is the number of timesteps in the data sample.
Ut=160;

%Cell-arrays to store measurement data. 
M0zt=cell(1,Ut);
M1zt=cell(1,Ut);
M0xt=cell(1,Ut);
M1xt=cell(1,Ut);
M0yt=cell(1,Ut);
M1yt=cell(1,Ut);

for ttt=1:Ut
%The data files for each time step is in the same form as for the SPAM
%extraction. Please refer to the description of the format in
%'SPAM_extract.m' for details.    
    
%Load qubit file
load(['../Data/Data2020_2/qubit1_0/t_' num2str((ttt)) '.mat']);

%The vectors M0x,M1x,.. contains the meassurement results in the x,z,y-basis for
%the input vectors. 1: |0>, 2: |1>, 3: |+>, 4: |->, 5: |+i>, 5: |-i>
M0x(1,1)=(length(m0x)-sum(m0x));
M0x(1,2)=(length(m1x)-sum(m1x));
M0x(1,3)=(length(mpx)-sum(mpx));
M0x(1,4)=(length(mmx)-sum(mmx));
M0x(1,5)=(length(mipx)-sum(mipx));
M0x(1,6)=(length(mimx)-sum(mimx));

M1x(1,1)=sum(m0x);
M1x(1,2)=sum(m1x);
M1x(1,3)=sum(mpx);
M1x(1,4)=sum(mmx);
M1x(1,5)=sum(mipx);
M1x(1,6)=sum(mimx);

M0y(1,1)=(length(m0y)-sum(m0y));
M0y(1,2)=(length(m1y)-sum(m1y));
M0y(1,3)=(length(mpy)-sum(mpy));
M0y(1,4)=(length(mmy)-sum(mmy));
M0y(1,5)=(length(mipy)-sum(mipy));
M0y(1,6)=(length(mimy)-sum(mimy));

M1y(1,1)=sum(m0y);
M1y(1,2)=sum(m1y);
M1y(1,3)=sum(mpy);
M1y(1,4)=sum(mmy);
M1y(1,5)=sum(mipy);
M1y(1,6)=sum(mimy);

M0z(1,1)=(length(m0z)-sum(m0z));
M0z(1,2)=(length(m1z)-sum(m1z));
M0z(1,3)=(length(mpz)-sum(mpz));
M0z(1,4)=(length(mmz)-sum(mmz));
M0z(1,5)=(length(mipz)-sum(mipz));
M0z(1,6)=(length(mimz)-sum(mimz));

M1z(1,1)=sum(m0z);
M1z(1,2)=sum(m1z);
M1z(1,3)=sum(mpz);
M1z(1,4)=sum(mmz);
M1z(1,5)=sum(mipz);
M1z(1,6)=sum(mimz);

M0zt{1,ttt}=M0z;
M1zt{1,ttt}=M1z;
M0xt{1,ttt}=M0x;
M1xt{1,ttt}=M1x;
M0yt{1,ttt}=M0y;
M1yt{1,ttt}=M1y;

end

%The singel qubit Pauli-martrices: 
S{1,1}=1/sqrt(2)*[1 0;0 -1];
S{1,2}=1/sqrt(2)*[0,1;1,0];
S{1,3}=1/sqrt(2)*[0,-1i;1i,0];

%The initial guess for hnm is in this case based on a dephasing and amplitude damping
%channel: 

T1=10000;          %T1 time microsec
T2=60;          %T2 time microsec

hnm0p=zeros(3,3);
hnm0p(1,1)=1/2*1/T2;
hnm0p(2,2)=1/4*1/T1;
hnm0p(3,3)=1/4*1/T1;
hnm0p(2,3)=-1/4*1i*1/T1;
hnm0p(3,2)=conj(hnm0p(2,3));

hnm0p=nearestSPD(hnm0p);

hnm0=ctranspose(chol(hnm0p+1e-13*eye(size(hnm0p))));

%Transforming the matrix hnm0 to a vector for fminsearch. 
odd=[1:2:31];
x0=zeros(1,13);
for ttt=1:3
for jjj=ttt:3
if ttt==jjj
x0(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1)=hnm0(jjj,ttt);
else
x0(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1+2*(jjj-ttt))=imag(hnm0(jjj,ttt));
x0(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+2*(jjj-ttt))=real(hnm0(jjj,ttt));
end
end
end

%The initial guess for the Hamiltonian
H0=[0 0;0 0];
%Put in the right format for fminsearch
x0(1,10)=H0(1,1);
x0(1,11)=real(H0(1,2));
x0(1,12)=imag(H0(1,2));
x0(1,13)=H0(2,2);

%fminsearch is used since there are no constraints for the optimization. 
options = optimset('MaxFunEvals',30000);
[xm,L0]=fminsearch(@(x)likelihood_Lindblad(M0zt,M1zt,M0xt,M1xt,M0yt,M1yt,POVM,S,rho,x),x0,options);
x0=xm;

%maximum log-likelihood
maxL=-L0;

%Extract the Hamiltonian
H(1,1)=xm(1,10);
H(1,2)=xm(1,11)+1i*xm(1,12);
H(2,1)=conj(H(1,2));
H(2,2)=xm(1,13);


%Extract the Lindblad matrix
hnmp=zeros(3,3);
for ttt=1:3
for jjj=ttt:3
if ttt==jjj
hnmp(jjj,ttt)=abs(xm(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1));
else
hnmp(jjj,ttt)=xm(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+2*(jjj-ttt))+1i*xm(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1+2*(jjj-ttt));
end
end
end

hnm=hnmp*ctranspose(hnmp);

%% Get the Lindblad operators and decay rates: 

%pauli operators
s1=[1,0;0,-1];
s2=[0,1;1,0];
s3=[0,-1i;1i,0];

[V,D]=eig(hnm);

L1=V(1,1)*s1+V(2,1)*s2+V(3,1)*s3;
L2=V(1,2)*s1+V(2,2)*s2+V(3,2)*s3;
L3=V(1,3)*s1+V(2,3)*s2+V(3,3)*s3;

%Renormalization of jump operators.
n1=trace(L1*ctranspose(L1));
n2=trace(L2*ctranspose(L2));
n3=trace(L3*ctranspose(L3));

D(1,1)=D(1,1)*n1;
D(2,2)=D(2,2)*n1;
D(3,3)=D(3,3)*n1;

L1=L1./sqrt(n1);
L2=L2./sqrt(n2);
L3=L3./sqrt(n3);

%Save output to file
save('Lindblad2020_1_0.mat','L1','L2','L3','D','H','hnm','maxL','-mat')

toc;