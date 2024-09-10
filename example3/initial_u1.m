function result=initial_u1(X,Y,T,alpha,rho)
rho=30;
Y=Y';
Nx=length(X); Ny=length(Y); 
result=zeros(Ny,Nx);
line=find(Y<=0.5); 
result(line,:)=kron(ones(1,Nx),tanh(rho*(Y(line)-0.25)));

line=find(Y>0.5);
result(line,:)=kron(ones(1,Nx),tanh(rho*(0.75-Y(line))));

    