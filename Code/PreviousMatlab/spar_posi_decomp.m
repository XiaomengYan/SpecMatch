% sparse positive decomposition algorithm
% using the function: starlet second-generation decompostion
% starlet second-generation reconstruction
% function: get_m
function alpha=spar_posi_decomp(X)
N=size(X);
c=sec_starlet_trans(X);
N=size(c,2);
J=size(c,1);
%M=get_m2(c);
M=ones(J,N);
for i=1:J
    for j=1:N
        if i<=3
    M(i,j)=0;
        end
    end
end
lambda0=max(c(:));
alpha=zeros(J,N);
alphay=zeros(J,N);
for i=1:1000
    X=posi_reconstr(alpha);
    alphay=sec_starlet_trans(X);
    r=M.*(c-alphay);
    lambda=lambda0*(1-i/1000);
    alpha=alpha+r-lambda*ones(J,N);
    for j=1:J
        for n=1:N
            if alpha(j,n)>=0
                alpha(j,n)=alpha(j,n);
            end
            if alpha(j,n)<0
                alpha(j,n)=0;
            end
        end
    end
    i=i+1;
end
    
    
