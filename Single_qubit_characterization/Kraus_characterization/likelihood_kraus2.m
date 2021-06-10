function [Out] = likelihood_kraus2(f1,f2,f3,f4,f5,f6,M,K1,K2,K3,K4,rho,c)

%the single qubit Gates for qubit preparation:
sigmax=[0 -1i;-1i 0];
hadp=1/sqrt(2)*[1 -1;1 1];
hadip=1/sqrt(2)*[1 1i; 1i 1];
hadm=1/sqrt(2)*[1 -1;1 1];
hadim=1/sqrt(2)*[1 -1i; -1i 1];


%Get the evolved density matrix from the Kraus operators
p0=evolution(K1,K2,K3,K4,rho);
p1=evolution(K1,K2,K3,K4,sigmax*rho*ctranspose(sigmax));
pp=evolution(K1,K2,K3,K4,hadp*rho*ctranspose(hadp));
pm=evolution(K1,K2,K3,K4,hadm*rho*ctranspose(hadm));
pip=evolution(K1,K2,K3,K4,hadip*rho*ctranspose(hadip));
pim=evolution(K1,K2,K3,K4,hadim*rho*ctranspose(hadim));

if c==1
%Lz
Mb=M;
    Out=f1*log(trace(Mb*p0))+f2*log(trace(Mb*p1))+f3*log(trace(Mb*pp))+f4*log(trace(Mb*pm))+f5*log(trace(Mb*pip))+f6*log(trace(Mb*pim));
end
if c==2
%Lx
Mx=1/sqrt(2)*[1 1;-1 1];
Mb=ctranspose(Mx)*M*Mx;
    Out=f1*log(trace(Mb*p0))+f2*log(trace(Mb*p1))+f3*log(trace(Mb*pp))+f4*log(trace(Mb*pm))+f5*log(trace(Mb*pip))+f6*log(trace(Mb*pim));
end
if c==3
%Ly
My=1/sqrt(2)*[1 -1i; -1i 1];
Mb=ctranspose(My)*M*My;
    Out=f1*log(trace(Mb*p0))+f2*log(trace(Mb*p1))+f3*log(trace(Mb*pp))+f4*log(trace(Mb*pm))+f5*log(trace(Mb*pip))+f6*log(trace(Mb*pim));
end

%Due to limited numerical precision matlab can give NaN or InF numbers
%which would kill the optimization. The fix here is to give a large penalty
%for such points but keep the optimization running. A warning is at the
%same time being issued. 
if isnan(Out)
Out=-sum(sum(f));
disp('NaN in log-likelihood');
end
if isinf(Out)
Out=-sum(sum(f));
disp('InF in log-likelihood');
end

end

