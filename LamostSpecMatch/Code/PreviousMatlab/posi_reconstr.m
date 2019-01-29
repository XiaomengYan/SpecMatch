% Starlet Second-Generation Reconstruction Algorithm
% Postive dictionay
function S=posi_reconstr(X)
N=size(X,2);
J=size(X,1);
S=zeros(1,N);
h_1D=1/16.*[1,4,6,4,1];
c=zeros(J+1,N);
c(J,:)=X(J,:);
v=zeros(1,N);
for l=J-1:-1:1
     for j=1:N
        sum=0;
        for k=-2:1:2
            jpk=j+2^(l-1)*k;
            while jpk<1|jpk>N
                if jpk<1
                    jpk=2-jpk;
                end
                if jpk>N
                    jpk=2*N-jpk;
                end
            end
                %fprintf('%d %d %d\n',l-1,jpk);
                e=c(l+1,jpk)*h_1D(k+3);
                sum=sum+e;
        end
            v(1,j)=sum;
     end
     c(l,:)=v+X(l,:);
end
S=c(1,:);
    
    

