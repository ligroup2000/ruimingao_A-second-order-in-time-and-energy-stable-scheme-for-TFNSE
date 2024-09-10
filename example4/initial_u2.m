function result=initial_u2(X,Y,Z,T,alpha)
Nx=length(X); Ny=length(Y); Nz=length(Z);
M=length(T);
result=zeros(Nx,Ny,Nz,M);
Xmesh=kron(ones(Ny,1),X); Xmesh=repmat(Xmesh,1,1,Nz);
Ymesh=kron(ones(1,Nx),Y'); Ymesh=repmat(Ymesh,1,1,Nz);
Zmesh=reshape(Z,1,1,Nz); Zmesh=repmat(Zmesh,Ny,Nx,1);
result=-cos(Xmesh).*sin(Ymesh).*cos(Zmesh);
end
    