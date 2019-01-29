% Pyramidal Median transform 
%input the original signal- continuum 
%output the coefficient matrix and a length vector which records the length
%of every row of coefficient matrix
function [c,v]=pry_medi_tans(X)
N=size(X,2);
J=floor(log2(N));
c=zeros(J+1,N);%coefficient matrix w=[w1...wJ,cJ]
c(1,:)=X;% the zero level
c_med=zeros(J+1,N);
c_intp=zeros(J+1,N);
x_test=zeros(1,7);
w=zeros(J,N);% detail coefficient matrix 
n=N; % original dimension of the matrix (the first row)
v=zeros(1,J+1);
v(1,1)=N;
for l=2:J+1
    for j=1:n
        % the median of window 7 
    for i=j-3:j+3
        h=i;
        while h<1|h>n
            if h<1
                h=2-h;
            end
            if h>n
                h=2*n-h;
            end
        end 
        x_test(i+4-j)=c(l-1,h);
    end
    c_med(l,j)=median(x_test);
    end
    % obtain the odd term of c_med(j+1,:) to get c(j+1,:)
    v(1,l)=floor(n/2)+mod(n,2);
    for k=1:v(1,l)
        c(l,k)=c_med(l,2*k-1);
    end
    %interpolation
    t=2*v(1,l)-mod(n,2);
            % justfy n is odd or even
            % if n is odd, interpolte as usual
            if mod(t,2)==1
                for m=1:t
                    for k=1:v(1,l)
                        if m==2*k-1
                            c_intp(l,m)=c(l,k);
                        end
                        if m==2*k&&m<t
                            c_intp(l,m)=median([c(l,k),c(l,k+1)]);
                        end
                    end
                end
            end
            % if n is even, let the last term= the term before the last
            % term
            if mod(t,2)==0
                for m=1:t
                    for k=1:v(1,l)
                        if m==2*k-1
                            c_intp(l,m)=c(l,k);
                        end
                        if m==2*k&&m<t
                            c_intp(l,m)=median([c(l,k),c(l,k+1)]);
                        end
                    end
                end
                c_intp(l,t)=c_intp(l,t-1);
            end
    n=v(1,l);
    w(l-1,:)=c(l-1,:)-c_intp(l,:);
end
for l=1:J+1
    if l<=J
    c(l,:)=w(l,:);
    end
    if l==J+1
        c(l,:)=c(l,:);
    end
end
