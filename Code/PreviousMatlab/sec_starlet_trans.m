%Starlet Second-Generation Transform Algorithm
function c=sec_starlet_trans(X)
N=size(X,2);
J=floor(log2(N));
c=zeros(J+1,N);
v=zeros(J+1,N);
c(1,:)=X;% the zero level
w=zeros(J,N);
% define fiter h
h_1D=1/16.*[1,4,6,4,1];
% 1-D starlet Transform 
for l=2:J+1 %level
    % to get the analysising function coefficients
    %do UWT (with hole) towars each row and each entry of the row
    for j=1:N
        sum=0;
        for k=-2:1:2
            jpk=j+2^(l-2)*k;
            while jpk<1|jpk>N
                if jpk<1
                    jpk=2-jpk;
                end
                if jpk>N
                    jpk=2*N-jpk;
                end
            end
                %fprintf('%d %d %d\n',l-1,jpk);
                e=c(l-1,jpk)*h_1D(k+3);
                sum=sum+e;
        end
            c(l,j)=sum;
    end
    % second decomposition towards c_j level 
    for j=1:N
        sum=0;
        for k=-2:1:2
            jpk=j+2^(l-2)*k;
            % periodic&&reflexove bound
            while jpk<1|jpk>N
                if jpk<1
                    jpk=2-jpk;
                end
                if jpk>N
                    jpk=2*N-jpk;
                end
            end
            %fprintf('%d %d %d\n',l-1,jpk);
            e=c(l,jpk)*h_1D(k+3);
            sum=sum+e;
        end
            v(l,j)=sum;
    end   
    w(l,:)=c(l-1,:)-v(l,:);
end
for l=1:J+1
    if l<=J
    c(l,:)=w(l+1,:);
    end
    if l==J+1
        c(l,:)=c(l,:);
    end
end