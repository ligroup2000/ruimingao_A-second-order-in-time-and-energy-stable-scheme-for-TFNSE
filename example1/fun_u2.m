function result=fun_u2(X,Y,T,alpha)
Nx=length(X); Ny=length(Y); M=length(T);
result=zeros(Nx,Ny,M);
Xmesh=kron(ones(Ny,1),X);
Ymesh=kron(ones(1,Nx),Y');
A=-(sin(2*Xmesh)).*(sin(Ymesh).^2);
for k=1:M
    result(:,:,k)=(power(T(k),2)+power(T(k),alpha)+1)*A;
end
    
end