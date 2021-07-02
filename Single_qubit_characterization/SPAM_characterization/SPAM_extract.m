%% Script to obtain the SPAM errors
%In this script we use a maximum likelihood to extract the SPAM errors for one qubit 

%Load qubit data file giving measurement results for a single qubit for t=0. (path is
%arbitrary)
load('Data/Data2020_2/qubit1_p/t_0.mat');

%In our case the measurement data was stored as arrays like m0x meaning input state
%|0> and measurement basis x. m0x would then be a vector of 0s and 1s
%recording all the measruement outcomes. 1=|1>, p=|+>, m=|->, ip=|+i>,
%im=|-i>. 


%The vectors M0x,M1x,.. should contain the meassurement results in the x,z,y-basis for
%the input vectors. 1: |0>, 2: |1>, 3: |+>, 4: |->, 5: |+i>, 5: |-i>. Thus, 
%M0x (M1x) is then the number of recorded 0s (1s) in m0x. 

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

%% Optimization

Sample=10000; %number of different starting points for optimization
XM=cell(1,Sample); %will store the output vector of optimization
L=zeros(1,Sample); %will store the output loglikelihood

%Fmincon specifications such as constraints (see constraintpovm file)
nonlcon=@constraintpovm;
options = optimoptions('fmincon','MaxFunctionEvaluations',20000);
for kk=1:Sample
%Sample a starting point for the POVM (using cholesky decomposition)
POVM0p=samplepovm(0.1);
POVM0=ctranspose(chol(POVM0p));

x0(1,1)=POVM0(1,1);
x0(1,2)=real(POVM0(2,1));
x0(1,3)=imag(POVM0(2,1));
x0(1,4)=POVM0(2,2);

%Sample starting point for initial state rho (using cholesky decomposition)
rho0=samplerho(0.1);
rho0p=ctranspose(chol(rho0+1e-13*eye(size(rho0))));

x0(1,5)=rho0p(1,1);
x0(1,6)=rho0p(2,2);
x0(1,7)=real(rho0p(2,1));
x0(1,8)=imag(rho0p(2,1));

%The maximum likelihood optimization using fmincon optimization
[xm,L0]=fmincon(@(X)likelihood_povm1(M0z,M1z,M0x,M1x,M0y,M1y,X),x0,[],[],[],[],[],[],nonlcon,options);
L(1,kk)=-L0;
XM{1,kk}=xm;
end
[maxL,Lt]=max(L); %find maximum loglikelihood.
nummaxL=find(L==maxL);
if length(nummaxL)>1
disp('multiple minima'); %To give information about multiple maximas found. 
end
xm=XM{1,Lt}; %Optimal solution

%Extracting POVM from solution
POVMo(1,1)=xm(1,1);
POVMo(2,1)=xm(1,2)+1i*xm(1,3);
POVMo(2,2)=xm(1,4);

M0=POVMo*ctranspose(POVMo);


POVM=M0;

%Extracting the initial state density matrix from solution:
rhoo(1,1)=xm(1,5);
rhoo(2,1)=xm(1,7)+1i*xm(1,8);
rhoo(2,2)=xm(1,6);

rhop=rhoo*ctranspose(rhoo);
rho=1/trace(rhop)*rhop;

%save output to file
save('SPAM2020_1_p.mat','POVM','rho','minL','-mat');