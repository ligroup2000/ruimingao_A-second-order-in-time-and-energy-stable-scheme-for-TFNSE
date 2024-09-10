function result=f3(t,X,Y,Z,alpha,nu)
Nx=length(X); Ny=length(Y); Nz=length(Z);
result=zeros(Nx,Ny,Nz);
Xmesh=kron(ones(Ny,1),X); Xmesh=repmat(Xmesh,1,1,Nz);
Ymesh=kron(ones(1,Nx),Y'); Ymesh=repmat(Ymesh,1,1,Nz);
Zmesh=reshape(Z,1,1,Nz); Zmesh=repmat(Zmesh,Ny,Nx,1);

result=-2*(2/gamma(3-alpha)*power(t,2-alpha)+gamma(1+alpha))*cos(Xmesh).*cos(Ymesh).*sin(Zmesh)...
    -6*nu*(t^alpha+t^2+1)*cos(Xmesh).*cos(Ymesh).*sin(Zmesh)+...
    +2*power(t^alpha+t^2+1,2)*sin(Zmesh).*cos(Zmesh).*...
    ((sin(Xmesh).^2).*(cos(Ymesh).^2)+(cos(Xmesh).^2).*(sin(Ymesh).^2)+2*(cos(Xmesh).^2).*(cos(Ymesh).^2))...
    +(t^alpha+t^2+1)*sin(Xmesh).*sin(Ymesh).*cos(Zmesh);

end