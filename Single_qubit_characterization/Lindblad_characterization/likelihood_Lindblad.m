function [Out] = likelihood_Lindblad(M0zt,M1zt,M0xt,M1xt,M0yt,M1yt,POVM,S,rho,x)
%Maximum likelihood functional for the Lindblad extraction. 

%number of time steps in data file
U=160;

%get the Hamiltonian
H(1,1)=x(1,10);
H(1,2)=x(1,11)+1i*x(1,12);
H(2,1)=conj(H(1,2));
H(2,2)=x(1,13);


%Get the Lindblad matrix
odd=[1:2:31];
hnmp=zeros(3,3);
for ttt=1:3
for jjj=ttt:3
if ttt==jjj
hnmp(jjj,ttt)=abs(x(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1));
else
hnmp(jjj,ttt)=x(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+2*(jjj-ttt))+1i*x(1,6*(ttt-1)-sum(odd(1,1:(ttt-1)))+1+2*(jjj-ttt));
end
end
end

hnm=hnmp*ctranspose(hnmp);


%The single qubit rotations used for state preparation. 
sigmax=[0 -1i;-1i 0];
hadp=1/sqrt(2)*[1 -1;1 1];
hadip=1/sqrt(2)*[1 1i; 1i 1];
hadm=1/sqrt(2)*[1 -1;1 1];
hadim=1/sqrt(2)*[1 -1i; -1i 1];


%The initial densitry matrices. 
y{1,1}=rho;
y{1,2}=sigmax*rho*ctranspose(sigmax);
y{1,3}=hadp*rho*ctranspose(hadp);
y{1,4}=hadm*rho*ctranspose(hadm);
y{1,5}=hadip*rho*ctranspose(hadip);
y{1,6}=hadim*rho*ctranspose(hadim);

%Load the data time vector and use unit of 1 microsec.
load('../Data/Data2020_2/time.mat');
time2=time(1,1:161)*10^(6); 
time=time2;

L=zeros(1,U);
%Construct liouvillian superoperator from the Lindblad matrix

Id=eye(2);
Lou=zeros(4,4);
for ii=1:3
for jj=1:3
Lp=hnm(ii,jj).*(kron(S{1,ii},Id)*kron(Id,transpose(ctranspose(S{1,jj})))-1/2*(kron(ctranspose(S{1,jj})*S{1,ii},Id)+kron(Id,transpose(ctranspose(S{1,jj})*S{1,ii}))));    
Lou=Lou+Lp;
end
end

% Add the Hamiltonian evolution
Lou=1i*(kron(Id,transpose(H))-kron(H,Id))+Lou;

%Get the density matrices after evolution under the Louvillian
rhop=cell(1,5);
for kk=1:6
Yp=zeros(length(time),4);
for nn=1:length(time)
Yp(nn,:)=expm(Lou*time(1,nn))*reshape(y{1,kk},4,1); 
end
rhop{1,kk}=Yp;
end


%Calculate the log-likelihood for each timestep. 
for ii=1:U   
for kk=1:6
Ypp=rhop{1,kk};    
y{ii+1,kk}=reshape(Ypp(ii+1,:),2,2);
end

M0z=M0zt{1,ii};
M1z=M1zt{1,ii};
M0x=M0xt{1,ii};
M1x=M1xt{1,ii};
M0y=M0yt{1,ii};
M1y=M1yt{1,ii};

p0f=y{ii+1,1};
p1f=y{ii+1,2};
ppf=y{ii+1,3};
pmf=y{ii+1,4};
pipf=y{ii+1,5};
pimf=y{ii+1,6};


%The log-likelihood

Lz=likelihood_Lindblad2(M0z,POVM,p0f,p1f,ppf,pmf,pipf,pimf,1)+likelihood_Lindblad2(M1z,eye(2,2)-POVM,p0f,p1f,ppf,pmf,pipf,pimf,1);

Lx=likelihood_Lindblad2(M0x,POVM,p0f,p1f,ppf,pmf,pipf,pimf,2)+likelihood_Lindblad2(M1x,eye(2,2)-POVM,p0f,p1f,ppf,pmf,pipf,pimf,2);

Ly=likelihood_Lindblad2(M0y,POVM,p0f,p1f,ppf,pmf,pipf,pimf,3)+likelihood_Lindblad2(M1y,eye(2,2)-POVM,p0f,p1f,ppf,pmf,pipf,pimf,3);

L(1,ii)=-(Lz+Lx+Ly);

%Due to numerical imprecision Matlab can give imaginary output. 
%The fix below is to give warning about this if the imaginary error is too
%large and stear the optimization away from such points. Can be adjusted. 
if abs(imag(L(1,ii)))<10^(-10) 
L(1,ii)=real(L(1,ii));
else
L(1,ii)=abs(L(1,ii)*100000);
disp('imaginary output');
end
end

%The total -(log-likelihood)
Out=sum(L);
end

