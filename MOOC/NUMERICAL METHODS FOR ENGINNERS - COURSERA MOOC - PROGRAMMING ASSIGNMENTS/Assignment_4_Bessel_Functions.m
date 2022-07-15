% Defining number of roots and number of fucntions whose root is to be
% calculated

num_roots=5; num_functions=6;
%initial guess for roots (from Wolfram MathWorld)
zeros_guess=[2.4,3.8,5.1,6,7.5,8.7;...
   5.5,7,8.4,9.7,11,12;...
   8.6 10,11.6,13,14,16;...
   11.8,13,15,16,18,19;...
   15,16.4,18,19.4,21,22]; 
%Computing first five roots of first six Bessel functions
%Putting in variable bzeros with size(bzeros) = [5, 6]

f = @(x,n,theta) (1/pi) * (cos(x.*sin(theta) - n*theta));

bzeros = zeros(num_roots,num_functions);

%Computing the roots%
for k=1:num_roots
    for n=1:num_functions
        j = @(x) integral(@(theta)f(x,n-1,theta),0,pi);
        r = fzero(@(x)j(x),zeros_guess(k,n));
        bzeros(k,n)=r;
    end
end



%print table
fprintf('k     J0(x)     J1(x)     J2(x)     J3(x)     J4(x)     J5(x)\n')
for k=1:num_roots
    fprintf('%i',k)
    for n=0:num_functions-1
        fprintf('%10.4f',bzeros(k,n+1));
    end
    fprintf('\n');
end
