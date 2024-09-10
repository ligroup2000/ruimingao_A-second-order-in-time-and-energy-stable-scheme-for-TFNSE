function result=fun_u3(X,Y,Z,T,alpha)
Nx=length(X); Ny=length(Y); Nz=length(Z);
M=length(T);
result=zeros(Nx,Ny,Nz,M);
Xmesh=kron(ones(Ny,1),X); Xmesh=repmat(Xmesh,1,1,Nz);
Ymesh=kron(ones(1,Nx),Y'); Ymesh=repmat(Ymesh,1,1,Nz);
Zmesh=reshape(Z,1,1,Nz); Zmesh=repmat(Zmesh,Ny,Nx,1);
A=cos(Xmesh).*cos(Ymesh).*sin(Zmesh);
for k=1:M
    result(:,:,:,k)=-2*(power(T(k),2)+power(T(k),alpha)+1)*A;
end

end
    