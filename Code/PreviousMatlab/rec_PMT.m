function X=rec_PMT(c,v)
N=size(c,2);
J=size(c,1);%12
X=zeros(1,N);
c_intp=zeros(J,N);
%M=get_m(c);
%c=M.*c;
for l=J:-1:2
    n=v(1,l)*2-mod(v(1,l-1),2);
    if mod(v(1,l-1),2)==1
    for j=1:n
        for k=1:v(1,l)
             % justfy n is odd or even
            % if n is odd, interpolte as usual
        if j==2*k-1
        c_intp(l,j)=c(l,k);
        end
        if j==2*k&&j<n
            c_intp(l,j)=median([c(l,k),c(l,k+1)]);
        end
        end
    end
    end
    % if n is even, let the last term= the term before the last
    % term
    if mod(v(1,l-1),2)==0
        for j=1:n
            for k=1:v(1,l)
                if j==2*k-1
                    c_intp(l,j)=c(l,k);
                end
                if j==2*k&&j<n
                    c_intp(l,j)=median([c(l,k),c(l,k+1)]);
                end
            end
        end
        c_intp(l,n)=c_intp(l,n-1);
    end
    c(l-1,:)=c(l-1,:)+c_intp(l,:);
end
X=c(1,:);

    