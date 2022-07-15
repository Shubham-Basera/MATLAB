clc;clear;

% Defining paramters
r=28; sigma=10; beta=8/3; 
x1=0; y1=0; z1=0;
x2=sqrt(beta*(r-1)); y2=sqrt(beta*(r-1)); z2=r-1;
x3=-sqrt(beta*(r-1)); y3=-sqrt(beta*(r-1)); z3=r-1;
nx=500; nz=500;
xmin=-40; xmax=40; zmin=-40; zmax=40;
x_grid=linspace(xmin,xmax,nx); z_grid=linspace(zmin,zmax,nz);
[X,Z]=meshgrid(x_grid,z_grid);

% WritING Newton's method using every grid point as the initial condition
RelTol=1.e-06; AbsTol=1.e-09;
for ix=1:nx
    for iz=1:nz
        x=X(iz,ix); y=3*sqrt(2); z=Z(iz,ix); % Setting initial value y=3*sqrt(2) for all values (x,z) on the grid
        error=Inf;
        while error > max(RelTol*max(abs([x,y,z])),AbsTol)
            J=[-sigma sigma 0;(r-z),-1,-x;y,x,-beta]; % DEFINING THE JACOBIAN MATRIX
            rhs = -[sigma*(y-x);r*x-y-x*z;x*y-beta*z]; % DEFINING THE RIGHT-HAND SIDE
            delta_xyz=J\rhs;
            x = x + delta_xyz(1);
            y = y + delta_xyz(2);
            z = z + delta_xyz(3);
            error=max(abs(delta_xyz));
        end
        X(iz,ix)=x; Z(iz,ix)=z; % % Saving the values of the converged roots

    end
end


% PLOTTING FRACTAL %
eps=1.e-03;
X1 = abs(X-x1) < eps; X2 = abs(X-x2) < eps; X3 = abs(X-x3) < eps;
X4 = ~(X1+X2+X3);
figure; 
map = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; colormap(map); %[red;green;blue;black]
X=(X1+2*X2+3*X3+4*X4); 
image([xmin xmax], [zmin zmax], X); set(gca,'YDir','normal');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize',14);
title('Fractal from the Lorenz Equations', 'Interpreter', 'latex','FontSize', 16)  
