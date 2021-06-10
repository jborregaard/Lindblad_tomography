%% Script to obtain the Kraus operators
%In this script we use a maximum likelihood to extract the Kraus operators 
tic;

%U is the number of timesteps in the data sample. 
U=160;

Kraus=cell(1,U); %Extracted Kraus operators will be stored here
maxL=zeros(1,U); %maximum log-likelihood will be stored here
for ii=1:U

%% Load data and reformat
%The data files for each time step is in the same form as for the SPAM
%extraction. Please refer to the description of the format in
%'SPAM_extract.m' for details.

%Load qubit file
load(['Data/Data2020_2/qubit1_1/t_' num2str((ii)) '.mat'])

%Load data time vector
load(['Data/Data2020_2/time.mat'])

%load extracted SPAM characterization
load('SPAM2020_1_1.mat');


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

%% Initial guess Krauss operators
if ii==1
%Initial guess of the Kraus operators i.e. for the first time step. For small initial time step this should be close to identity channel.  
K1a=[1 0;0 1];
K2a=[0 0;0 0];
K1d=[1 0;0 1];
K2d=[0 0;0 0];

Kp1=K1d*K1a;
Kp2=K1d*K2a;
Kp3=K2d*K1a;
Kp4=K2d*K2a;

else
%In future timesteps the Krauss operators from the previous step is used as initial guess.     
Kp1=Kraus{1,ii-1};
Kp2=Kraus{2,ii-1};
Kp3=Kraus{3,ii-1};
Kp4=Kraus{4,ii-1};  

end

%Bring Kraus operators in to vector for optimization
x0(1,1)=real(Kp1(1,1));
x0(1,2)=imag(Kp1(1,1));
x0(1,3)=real(Kp1(1,2));
x0(1,4)=imag(Kp1(1,2));
x0(1,5)=real(Kp1(2,1));
x0(1,6)=imag(Kp1(2,1));
x0(1,7)=real(Kp1(2,2));
x0(1,8)=imag(Kp1(2,2));

x0(1,9)=real(Kp2(1,1));
x0(1,10)=imag(Kp2(1,1));
x0(1,11)=real(Kp2(1,2));
x0(1,12)=imag(Kp2(1,2));
x0(1,13)=real(Kp2(2,1));
x0(1,14)=imag(Kp2(2,1));
x0(1,15)=real(Kp2(2,2));
x0(1,16)=imag(Kp2(2,2));

x0(1,17)=real(Kp3(1,1));
x0(1,18)=imag(Kp3(1,1));
x0(1,19)=real(Kp3(1,2));
x0(1,20)=imag(Kp3(1,2));
x0(1,21)=real(Kp3(2,1));
x0(1,22)=imag(Kp3(2,1));
x0(1,23)=real(Kp3(2,2));
x0(1,24)=imag(Kp3(2,2));

x0(1,25)=real(Kp4(1,1));
x0(1,26)=imag(Kp4(1,1));
x0(1,27)=real(Kp4(1,2));
x0(1,28)=imag(Kp4(1,2));
x0(1,29)=real(Kp4(2,1));
x0(1,30)=imag(Kp4(2,1));
x0(1,31)=real(Kp4(2,2));
x0(1,32)=imag(Kp4(2,2));



%% The optimization
%The log-likelihood that we wish to maximize is defined in the function
%likelihood_kraus1. We use fmincon for the optimization.

nonlcon=@constraintkraus; %contraint on Kraus operators
options = optimoptions('fmincon','MaxFunctionEvaluations',10000);
warning('off')
[xm,L0]=fmincon(@(X)likelihood_kraus1(M0z,M1z,M0x,M1x,M0y,M1y,POVM,rho,X),x0,[],[],[],[],[],[],nonlcon,options);

%maximum log-likelihood
maxL(1,ii)=-L0;

%The estimated Kraus operators
K1=[xm(1,1)+1i*xm(1,2) xm(1,3)+1i*xm(1,4);
    xm(1,5)+1i*xm(1,6) xm(1,7)+1i*xm(1,8)];

K2=[xm(1,9)+1i*xm(1,10) xm(1,11)+1i*xm(1,12);
    xm(1,13)+1i*xm(1,14) xm(1,15)+1i*xm(1,16)];

K3=[xm(1,17)+1i*xm(1,18) xm(1,19)+1i*xm(1,20);
    xm(1,21)+1i*xm(1,22) xm(1,23)+1i*xm(1,24)];

K4=[xm(1,25)+1i*xm(1,26) xm(1,27)+1i*xm(1,28);
    xm(1,29)+1i*xm(1,30) xm(1,31)+1i*xm(1,32)];


Kraus{1,ii}=K1;
Kraus{2,ii}=K2;
Kraus{3,ii}=K3;
Kraus{4,ii}=K4;

end

%Save output Krauss operators
save('Kraus2020_2_1_1.mat','Kraus','maxL','-mat');

toc;