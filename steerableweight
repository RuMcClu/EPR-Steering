%calculates steerable weight
rhoini=load('filename.mat');
%input matrix of 4x4 matrices in .mat format
rhofield=fieldnames(rhoini);
answer=zeros(1,length(rhofield));
S=6;
RHO=zeros(2,2,S);
Eig=zeros(2,2,S);
PROJ=zeros(2,2,S);
mixpar=1;

Dl=zeros(8,S);
Lambda=zeros(8,3);
%SigmaL=zeros(2,2,8);
%sigmaUS=zeros(2,2,6);
%initial state


x=[0 1;1 0];
y=-[0 1i;-1i 0];
z=[1 0; 0 -1];
e=[1 0; 0 1];
% Pauli matrices

[XX,Xnoneed]=eig(x);
[YY,Ynoneed]=eig(y);
[ZZ,Znoneed]=eig(z);

% Pauli eigenvectors

EigVec(:,:,1)=XX(1:end,1);
EigVec(:,:,2)=XX(1:end,2);
EigVec(:,:,3)=YY(1:end,1);
EigVec(:,:,4)=YY(1:end,2);
EigVec(:,:,5)=ZZ(1:end,1);
EigVec(:,:,6)=ZZ(1:end,2);
%Pauli projectors
EigVal(1)=Xnoneed(1,1);
EigVal(2)=Xnoneed(2,2);
EigVal(3)=Ynoneed(1,1);
EigVal(4)=Ynoneed(2,2);
EigVal(5)=Znoneed(1,1);
EigVal(6)=Znoneed(2,2);
dimen=[2,2];
%subsystem dimensions for TrX

MM=0;
for i=1:2,
    for j=1:2,
        for k=1:2,
            MM=MM+1;
            Lambda(MM,1)=round(EigVec(:,:,i)'*x*EigVec(:,:,i));
            Lambda(MM,2)=round(EigVec(:,:,j+2)'*y*EigVec(:,:,j+2));
            Lambda(MM,3)=EigVec(:,:,k+4)'*z*EigVec(:,:,k+4);
%classical random variables
        end
    end
end

for i=1:6,
    if  i==1||i==2
        Kk=1;
    elseif i==3||i==4
        Kk=2;
    else
        Kk=3;
    end


    for j=1:8,
        if  Lambda(j,Kk)==EigVal(i)
            Dl(j,i)=1;
        else
            Dl(j,i)=0;
        end
    end
end
%extremal deterministic single-party conditional probability distributions

for nnn=1:length(rhofield),
    rho=rhoini.(rhofield{nnn});
    for i=1:S,
        PROJ(:,:,i)=kron(EigVec(:,:,i),EigVec(:,:,i)');
        RHO(:,:,i)=TrX(kron(PROJ(:,:,i),eye(2))*rho,1,dimen);
%evolved states (sub-normalised)
    end

    cvx_begin quiet sdp
   %requires cvx package to work
    cvx_precision best
    variable Sigma1(2,2) hermitian semidefinite;
    variable Sigma2(2,2) hermitian semidefinite;
    variable Sigma3(2,2) hermitian semidefinite;
    variable Sigma4(2,2) hermitian semidefinite;
    variable Sigma5(2,2) hermitian semidefinite;
    variable Sigma6(2,2) hermitian semidefinite;
    variable Sigma7(2,2) hermitian semidefinite;
    variable Sigma8(2,2) hermitian semidefinite;
    maximize(trace(Sigma1+Sigma2+Sigma3+Sigma4+Sigma5+Sigma6+Sigma7+Sigma8));
    for k=1:6
        (RHO(:,:,k)-(Dl(1,k)*Sigma1+Dl(2,k)*Sigma2+Dl(3,k)*Sigma3+Dl(4,k)*Sigma4+Dl(5,k)*Sigma5+Dl(6,k)*Sigma6+Dl(7,k)*Sigma7+Dl(8,k)*Sigma8))==hermitian_semidefinite(2);
        %==semidefinite(2);
    end
    %for k=1:8
    %(SigmaL(:,:,k))>=0;
    %==semidefinite(2);
    %end
    trace(Sigma1+Sigma2+Sigma3+Sigma4+Sigma5+Sigma6+Sigma7+Sigma8)>=0;
    trace(Sigma1+Sigma2+Sigma3+Sigma4+Sigma5+Sigma6+Sigma7+Sigma8)<=1;
    cvx_end


    answer(nnn)=1-trace(Sigma1+Sigma2+Sigma3+Sigma4+Sigma5+Sigma6+Sigma7+Sigma8);
end
%export()
disp('finished')
