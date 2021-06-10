%% M-file to fit the Lindblad operators
%In this script we use a maximum likelihood to extract the Lindblad operators 

%Load extracted SPAM
load('../SPAM2020_2.mat');

%Ut is the number of timesteps in the data sample.
Ut=160;


%Cell-arrays to store measurement data. 
m00t=cell(1,Ut);
m01t=cell(1,Ut);
m10t=cell(1,Ut);
m11t=cell(1,Ut);
for ttt=1:Ut
    
%Load qubit file
%The data files for each time step is in the same form as for the SPAM
%extraction. Please refer to the description of the format in
%'SPAM_extract.m' for details.
load(['../Data/Data2020_2/Qb2/t_' num2str((ttt)) '.mat']);
load(['../Data/Data2020_2/Qb1/t_' num2str((ttt)) '.mat']);

U=1000; %Number of measurements per setting

%The measurement results are combined into # of 00, 01, 10, and 11 outcomes. (see SPAM_extract.m)  
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
m00t{1,ttt}=m00;
m01t{1,ttt}=m01;
m10t{1,ttt}=m10;
m11t{1,ttt}=m11;
end

%The 2-qubit operator-basis: 
ss{1,1}=1/sqrt(2)*[1 0;0 1];
ss{1,2}=1/sqrt(2)*[1 0;0 -1];
ss{1,3}=1/sqrt(2)*[0,1;1,0];
ss{1,4}=1/sqrt(2)*[0,-1i;1i,0];

S=cell(1,15);
for ii=1:4
for jj=1:4    
if jj==1 && ii==1
else    
S{1,4*(jj-1)+ii-1}=kron(ss{1,ii},ss{1,jj});
end
end
end


%As Initial guess for the Lindblad matrix, we can use the single qubit Lindblad matrices
%extracted for the qubits individually.

load('singlequbitLindblad.mat');

hnm0p=zeros(15,15);
for ii=1:3
for jj=1:3
if ii==jj
hnm1(ii,jj)=real(hnm1(ii,jj));
hnm2(ii,jj)=real(hnm2(ii,jj));
end
hnm0p(ii,jj)=hnm1(ii,jj);
hnm0p(4*ii,4*jj)=hnm2(ii,jj);
end
end



hnm0p=hnm;

%Bring Lindblad matricx into optim vector form (also using Cholesky
%decomposition). 
hnm0=ctranspose(chol(hnm0p+1e-13*eye(size(hnm0p))));

odd=[1:2:31];
x0=zeros(1,241);
for ttt=1:15
for jjj=ttt:15
if ttt==jjj
x0(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+1)=hnm0(jjj,ttt);
else
x0(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+1+2*(jjj-ttt))=imag(hnm0(jjj,ttt));
x0(1,30*(ttt-1)-sum(odd(1,1:(ttt-1)))+2*(jjj-ttt))=real(hnm0(jjj,ttt));
end
end
end

% Initial guess for Hamiltonian from extracted single qubit hamiltonians
H=kron(H1,H2);

% Input to optim vector
x0(1,226)=H(1,1);
x0(1,227)=real(H(1,2));
x0(1,228)=imag(H(1,2));
x0(1,229)=real(H(1,3));
x0(1,230)=imag(H(1,3));
x0(1,231)=real(H(1,4));
x0(1,232)=imag(H(1,4));
x0(1,233)=H(2,2);
x0(1,234)=real(H(2,3));
x0(1,235)=imag(H(2,3));
x0(1,236)=real(H(2,4));
x0(1,237)=imag(H(2,4));
x0(1,238)=H(3,3);
x0(1,239)=real(H(3,4));
x0(1,240)=imag(H(3,4));
x0(1,241)=H(4,4);


x0=real(x0);

%fminsearch is used since there are no constraints for the optimization. 
options = optimset('MaxFunEvals',6000); 
[xm,L0]=fminsearch(@(x)likelihood_Lindblad1(m00t,m01t,m10t,m11t,M00,M01,M10,M11,S,rho,x),x0,options);

%maximum log-likelihood
maxL=-L0;

% Get the extracted Hamiltonian
H(1,1)=xm(1,226);
H(1,2)=xm(1,227)+1i*xm(1,228);
H(2,1)=conj(H(1,2));
H(1,3)=xm(1,229)+1i*xm(1,230);
H(3,1)=conj(H(1,3));
H(1,4)=xm(1,231)+1i*xm(1,232);
H(4,1)=conj(H(1,4));
H(2,2)=xm(1,233);
H(2,3)=xm(1,234)+1i*xm(1,235);
H(3,2)=conj(H(2,3));
H(2,4)=xm(1,236)+1i*xm(1,237);
H(4,2)=conj(H(2,4));
H(3,3)=xm(1,238);
H(3,4)=xm(1,239)+1i*xm(1,240);
H(4,3)=conj(H(3,4));
H(4,4)=xm(1,241);

%Get the extracted Lindblad matrix
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

%Get the Lindblad operators and decay rates: 

[V,D]=eig(hnm);

L=cell(1,15);
for ii=1:15
L{1,ii}=V(1,ii)*S{1,1}+V(2,ii)*S{1,2}+V(3,ii)*S{1,3}+V(4,ii)*S{1,4}+V(5,ii)*S{1,5}+...
V(6,ii)*S{1,6}+V(7,ii)*S{1,7}+V(8,ii)*S{1,8}+V(9,ii)*S{1,9}+V(10,ii)*S{1,10}+...
V(11,ii)*S{1,11}+V(12,ii)*S{1,12}+V(13,ii)*S{1,13}+V(14,ii)*S{1,14}+V(15,ii)*S{1,15};
end


%Renormalization of jump operators and decay rates
for ii=1:15
n=trace(L{1,ii}*ctranspose(L{1,ii}));
D(1,ii)=D(1,ii)*n;
L{1,ii}=L{1,ii}/sqrt(n);
end

%Save output
save('Lindblad2020.mat','L','D','hnm','H', 'maxL','-mat')
