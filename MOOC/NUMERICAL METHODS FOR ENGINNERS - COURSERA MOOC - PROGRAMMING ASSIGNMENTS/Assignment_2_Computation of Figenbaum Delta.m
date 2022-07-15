clc;clear;

% Defining parameters
num_doublings=11; delta=zeros(1,num_doublings); delta(1)=5;
% Computing the Feigenbaum delta
m1=2; m2=1+sqrt(5);
for n = 2:num_doublings
    guess= m2 + (m2-m1)/delta(n-1);
    for i = 1:10
        x=1/2; xder=0;
        for j =1:2^(n)
            xder =x*(1-x) + guess*xder*(1-2*x);
            x =guess*x*(1-x);
        end
        guess = guess - (x-1/2)/xder;
    end
    m3= guess;
    delta(n)= (m2-m1)/(m3-m2);
    m1=m2; m2=m3;
end


% Outputing results
fprintf('n        delta(n)\n');
for n=1:num_doublings
    fprintf('%2g %18.15f\n',n,delta(n));
end
