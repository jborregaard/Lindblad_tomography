%% Script to obtain the Kraus operators
%In this script we use a maximum likelihood to extract the Kraus operators
%for 2-qubits

%Load extracted SPAM errors
load('SPAM2020_2.mat');

%Ut is the number of timesteps in the data sample.
Ut=160;

Kraus=cell(16,Ut);%Extracted Kraus operators will be stored here
maxL=zeros(1,Ut); %maximum log-likelihood will be stored here
for ttt=1:Ut 
    
%Load qubit file
%The data files for each time step is in the same form as for the SPAM
%extraction. Please refer to the description of the format in
%'SPAM_extract.m' for details.
load(['Data/Data2020_2/Qb2/t_' num2str((ttt)) '.mat']);
load(['Data/Data2020_2/Qb1/t_' num2str((ttt)) '.mat']);
U=1000; %Number of measurements per setting


%The measurement results are combined into # of 00, 01, 10, and 11 outcomes.(see SPAM_extract.m)  

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

if ttt==1
%Initial guess of the Kraus operators i.e. for the first time step. For small initial time step this should be close to identity channel.  
%Single qubit Kraus
K1a=[1 0;0 1];
K2a=[0 0;0 0];
K1d=[1 0;0 1];
K2d=[0 0;0 0];

Kp1=K1d*K1a;
Kp2=K1d*K2a;
Kp3=K2d*K1a;
Kp4=K2d*K2a;

%Two qubit Kraus
K11=kron(K1a,eye(2,2));
K12=kron(K2a,eye(2,2));
K13=kron(K1d,eye(2,2));
K14=kron(K2d,eye(2,2));

%Qubit 2
K21=kron(eye(2,2),K1a);
K22=kron(eye(2,2),K2a);
K23=kron(eye(2,2),K1d);
K24=kron(eye(2,2),K2d);

K01=K11*K21;
K02=K11*K22;
K03=K11*K23;
K04=K11*K24;
K05=K12*K21;
K06=K12*K22;
K07=K12*K23;
K08=K12*K24;
K09=K13*K21;
K010=K13*K22;
K011=K13*K23;
K012=K13*K24;
K013=K14*K21;
K014=K14*K22;
K015=K14*K23;
K016=K14*K24;

else
%In future timesteps the Krauss operators from the previous step is used as initial guess (Can be changed to output from previous optimizations)    

K01=Kraus{1,ttt-1};
K02=Kraus{2,ttt-1};
K03=Kraus{3,ttt-1};
K04=Kraus{4,ttt-1};
K05=Kraus{5,ttt-1};
K06=Kraus{6,ttt-1};
K07=Kraus{7,ttt-1};
K08=Kraus{8,ttt-1};
K09=Kraus{9,ttt-1};
K010=Kraus{10,ttt-1};
K011=Kraus{11,ttt-1};
K012=Kraus{12,ttt-1};
K013=Kraus{13,ttt-1};
K014=Kraus{14,ttt-1};
K015=Kraus{15,ttt-1};
K016=Kraus{16,ttt-1};
end

%Put Kraus operators into form of optimization vector
X0(1,1:16)=reshape(real(K01),1,16);
X0(1,1+16:2*16)=reshape(real(K02),1,16);
X0(1,1+2*16:3*16)=reshape(real(K03),1,16);
X0(1,1+3*16:4*16)=reshape(real(K04),1,16);
X0(1,1+4*16:5*16)=reshape(real(K05),1,16);
X0(1,1+5*16:6*16)=reshape(real(K06),1,16);
X0(1,1+6*16:7*16)=reshape(real(K07),1,16);
X0(1,1+7*16:8*16)=reshape(real(K08),1,16);
X0(1,1+8*16:9*16)=reshape(real(K09),1,16);
X0(1,1+9*16:10*16)=reshape(real(K010),1,16);
X0(1,1+10*16:11*16)=reshape(real(K011),1,16);
X0(1,1+11*16:12*16)=reshape(real(K012),1,16);
X0(1,1+12*16:13*16)=reshape(real(K013),1,16);
X0(1,1+13*16:14*16)=reshape(real(K014),1,16);
X0(1,1+14*16:15*16)=reshape(real(K015),1,16);
X0(1,1+15*16:16*16)=reshape(real(K016),1,16);

X0(1,256+1:256+16)=reshape(imag(K01),1,16);
X0(1,256+1+16:256+2*16)=reshape(imag(K02),1,16);
X0(1,256+1+2*16:256+3*16)=reshape(imag(K03),1,16);
X0(1,256+1+3*16:256+4*16)=reshape(imag(K04),1,16);
X0(1,256+1+4*16:256+5*16)=reshape(imag(K05),1,16);
X0(1,256+1+5*16:256+6*16)=reshape(imag(K06),1,16);
X0(1,256+1+6*16:256+7*16)=reshape(imag(K07),1,16);
X0(1,256+1+7*16:256+8*16)=reshape(imag(K08),1,16);
X0(1,256+1+8*16:256+9*16)=reshape(imag(K09),1,16);
X0(1,256+1+9*16:256+10*16)=reshape(imag(K010),1,16);
X0(1,256+1+10*16:256+11*16)=reshape(imag(K011),1,16);
X0(1,256+1+11*16:256+12*16)=reshape(imag(K012),1,16);
X0(1,256+1+12*16:256+13*16)=reshape(imag(K013),1,16);
X0(1,256+1+13*16:256+14*16)=reshape(imag(K014),1,16);
X0(1,256+1+14*16:256+15*16)=reshape(imag(K015),1,16);
X0(1,256+1+15*16:256+16*16)=reshape(imag(K016),1,16);

%% The optimization

%The log-likelihood that we wish to maximize is defined in the function
%likelihood_kraus1. We use fmincon for the optimization.
nonlcon=@constraint_kraus1; %Constraint on Kraus operators 
options = optimoptions('fmincon','MaxFunctionEvaluations',60000);

[Xm,L0]=fmincon(@(X)likelihood_Kraus1(m00,m01,m10,m11,M00,M01,M10,M11,rho,X),X0,[],[],[],[],[],[],nonlcon,options);

%maximum log-likelihood
maxL(1,ttt)=-L0;

% The estimated Kraus opeators
K1=reshape(Xm(1,1:16),4,4)+1i*reshape(Xm(1,256+1:256+16),4,4);
K2=reshape(Xm(1,1+16:2*16),4,4)+1i*reshape(Xm(1,256+1+16:256+2*16),4,4);
K3=reshape(Xm(1,1+2*16:3*16),4,4)+1i*reshape(Xm(1,256+1+2*16:256+3*16),4,4);
K4=reshape(Xm(1,1+3*16:4*16),4,4)+1i*reshape(Xm(1,256+1+3*16:256+4*16),4,4);
K5=reshape(Xm(1,1+4*16:5*16),4,4)+1i*reshape(Xm(1,256+1+4*16:256+5*16),4,4);
K6=reshape(Xm(1,1+5*16:6*16),4,4)+1i*reshape(Xm(1,256+1+5*16:256+6*16),4,4);
K7=reshape(Xm(1,1+6*16:7*16),4,4)+1i*reshape(Xm(1,256+1+6*16:256+7*16),4,4);
K8=reshape(Xm(1,1+7*16:8*16),4,4)+1i*reshape(Xm(1,256+1+7*16:256+8*16),4,4);
K9=reshape(Xm(1,1+8*16:9*16),4,4)+1i*reshape(Xm(1,256+1+8*16:256+9*16),4,4);
K10=reshape(Xm(1,1+9*16:10*16),4,4)+1i*reshape(Xm(1,256+1+9*16:256+10*16),4,4);
K11=reshape(Xm(1,1+10*16:11*16),4,4)+1i*reshape(Xm(1,256+1+10*16:256+11*16),4,4);
K12=reshape(Xm(1,1+11*16:12*16),4,4)+1i*reshape(Xm(1,256+1+11*16:256+12*16),4,4);
K13=reshape(Xm(1,1+12*16:13*16),4,4)+1i*reshape(Xm(1,256+1+12*16:256+13*16),4,4);
K14=reshape(Xm(1,1+13*16:14*16),4,4)+1i*reshape(Xm(1,256+1+13*16:256+14*16),4,4);
K15=reshape(Xm(1,1+14*16:15*16),4,4)+1i*reshape(Xm(1,256+1+14*16:256+15*16),4,4);
K16=reshape(Xm(1,1+15*16:16*16),4,4)+1i*reshape(Xm(1,256+1+15*16:256+16*16),4,4);

Kraus{1,ttt}=K1;
Kraus{2,ttt}=K2;
Kraus{3,ttt}=K3;
Kraus{4,ttt}=K4;
Kraus{5,ttt}=K5;
Kraus{6,ttt}=K6;
Kraus{7,ttt}=K7;
Kraus{8,ttt}=K8;
Kraus{9,ttt}=K9;
Kraus{10,ttt}=K10;
Kraus{11,ttt}=K11;
Kraus{12,ttt}=K12;
Kraus{13,ttt}=K13;
Kraus{14,ttt}=K14;
Kraus{15,ttt}=K15;
Kraus{16,ttt}=K16;


end

%Save output
save('Kraus2020_2.mat','Kraus','maxL','-mat');
