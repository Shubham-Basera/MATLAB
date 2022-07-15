clc; clear;

%%%%% Defining the square and grid parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=1;  %square is 2L x 2L 
N=100; %# of intervals in x and y directions
n=N+1; %# of gridpoints in x,y directions including boundaries
h=2*L/N;  %grid size in x,y directions
x=-L + (0:N)*h; %x values on the grid
y=-L + (0:N)*h; %y values on the grid
[X,Y]=meshgrid(x,y);

%%%%% Defining the indices associated with the boundaries %%%%%%%%%%%%%%%%%%%
% boundary_index = [bottom, left, top, right]
boundary_index=[          1:n,       1:n:1+(n-1)*n, ...
                1+(n-1)*n:n*n,   n:n:n*n           ]; 
%%%%% Diffusion constant and time-step parameters
D=1;
dt=h^2/(2*D); %borderline stability of FTCS scheme
alpha=dt*D/h^2; %equation parameter
nsteps=1000; %number of time steps
%%%%% CONSTRUCTing THE MATRIX AND COMPUTING LU DECOMPOSITION %%%%%%%%%%%%%%%%%%%%
diagonals = [ 2*(1+2*alpha)*ones(n*n,1),-alpha*ones(n*n,4)];
A = spdiags(diagonals, [0 -1 1 -n n], n*n, n*n); % use sparse matrices
I = speye(n*n);
A(boundary_index,:) = I(boundary_index,:);

% LU Decomposition%
[Low,U] = lu(A);

%%%%% Defining initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(n,n,nsteps);
sigma=L/4;
u(:,:,1)=(1/(2*pi*sigma^2)).*exp(-0.5*(X.^2+Y.^2)/sigma^2); 
u(1,:,1)=0; u(n,:,1)=0; u(:,1,1)=0; u(:,n,1)=0; %b.c.
%%%%% ADVANCING SOLUTION u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=2:nsteps
    b=zeros(n,n);
    b(2:n-1,2:n-1)=alpha*u(1:n-2,2:n-1,m-1) + 2*(1-2*alpha)*u(2:n-1,2:n-1,m-1) + alpha*u(3:n,2:n-1,m-1) + alpha*u(2:n-1,1:n-2,m-1) + alpha*u(2:n-1,3:n,m-1);
    b=reshape(b,n*n,1);       %column vector using natural ordering-same as b=b(:) 
    u(:,:,m)=reshape(U\(Low\b),n,n);
end
% %%%%% Plotting with animation%%%
figure('units','normalized','outerposition',[0 0 1 1])
s=surf(X,Y,u(:,:,1)); zlim([0, 2.6]);
xlabel('$x$','Interpreter','latex','FontSize',14); 
ylabel('$y$','Interpreter','latex','FontSize',14); 
zlabel('$u(x,y,t)$','Interpreter','latex','FontSize',14); 
title('Solution of the 2D diffusion equation','Interpreter','latex','FontSize',16);
pause(1)
for j=2:nsteps
     s.ZData=u(:,:,j); pause(0.01);
end
