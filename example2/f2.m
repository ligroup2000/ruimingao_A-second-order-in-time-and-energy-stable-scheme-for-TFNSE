function result=f2(t,X,Y,Z,alpha,nu)
Nx=length(X); Ny=length(Y); Nz=length(Z);
result=zeros(Nx,Ny,Nz);
Xmesh=kron(ones(Ny,1),X); Xmesh=repmat(Xmesh,1,1,Nz);
Ymesh=kron(ones(1,Nx),Y'); Ymesh=repmat(Ymesh,1,1,Nz);
Zmesh=reshape(Z,1,1,Nz); Zmesh=repmat(Zmesh,Ny,Nx,1);

result=(2/gamma(3-alpha)*power(t,2-alpha)+gamma(1+alpha))*cos(Xmesh).*sin(Ymesh).*cos(Zmesh)...
    +3*nu*(t^alpha+t^2+1)*cos(Xmesh).*sin(Ymesh).*cos(Zmesh)+...
    +power(t^alpha+t^2+1,2)*sin(Ymesh).*cos(Ymesh).*...
    ((cos(Xmesh).^2).*(cos(Zmesh).^2)-(sin(Xmesh).^2).*(cos(Zmesh).^2)+2*(cos(Xmesh).^2).*(sin(Zmesh).^2))...
    +(t^alpha+t^2+1)*sin(Xmesh).*cos(Ymesh).*sin(Zmesh);

end