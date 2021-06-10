function [Out] = likelihood_Lindblad2(f,M,p0,p1,pp,pm,pip,pim,c)
%Function with the lindblad log-likelihood functions
f11=f(1,1); f12=f(1,2); f13=f(1,3); f14=f(1,4); f15=f(1,5); f16=f(1,6);

if c==1
%Lz
Mb=M;
    Out=f11*log(trace(Mb*p0))+f12*log(trace(Mb*p1))+f13*log(trace(Mb*pp))+f14*log(trace(Mb*pm))+f15*log(trace(Mb*pip))+f16*log(trace(Mb*pim));
end
if c==2
%Lx
Ry=1/sqrt(2)*[1 1; -1 1];
Mb=ctranspose(Ry)*M*Ry;
    Out=f11*log(trace(Mb*p0))+f12*log(trace(Mb*p1))+f13*log(trace(Mb*pp))+f14*log(trace(Mb*pm))+f15*log(trace(Mb*pip))+f16*log(trace(Mb*pim));
end
if c==3
%Ly
Rx=1/sqrt(2)*[1 -1i; -1i 1];
Mb=ctranspose(Rx)*M*Rx;
    Out=f11*log(trace(Mb*p0))+f12*log(trace(Mb*p1))+f13*log(trace(Mb*pp))+f14*log(trace(Mb*pm))+f15*log(trace(Mb*pip))+f16*log(trace(Mb*pim));
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

