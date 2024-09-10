function result=initial_u2(X,Y,T,alpha)
delta=0.05;

Nx=length(X); Ny=length(Y); 
result=zeros(Ny,Nx);
result=kron(ones(Ny,1),delta*sin(2*pi*X));
    