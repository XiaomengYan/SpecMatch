%find the matrix corresponding the the multiresolution support
function M=get_m(wt)
% var(c0) and var(w1).
N=size(wt,2);
J=size(wt,1);
M=zeros(J,N);
wt_1=wt(2,:);
x_test=zeros(1,400);
sigma=zeros(J,N);
for l=1:N
    for i=l-200:l+200
        h=i;
        while h<1|h>N
            if h<1
                h=2-h;
            end
            if h>N
               h=2*N-h;
            end
        end
      x_test(i+201-l)=wt_1(h);
    end
        sigma(1,l)=median(abs(x_test-median(x_test)))/0.6745;
end
sigma0=sigma(1,:)*256/134;
% simulation to generate the variance of wl,l=2....J.
cov=zeros(N,N);
for i=1:N
    cov(i,i)=sigma0(i);
end
mu=wt(1,:);
p=mvnrnd(mu,cov,1000);
for i=1:1000
    B{i}=starlet_trans(p(i,:));
end
x=zeros(1,1000);
for j=2:J%w2:3.....
    for k=1:N
    for i=1:1000
        x(i)=B{i}(j,k);
    end
    sigma(j,k)=median(abs(x-median(x)))/0.6745;
    end
end
for i=1:6
    for j=1:N
        if abs(wt(i,j))/sigma(i,j)>4% sqrt(2logN)
            M(i,j)=1;
        end
    end
end
for i=7:J
    for j=1:N
        M(i,j)=1;
    end
end
    
    
    
    
    
    


            

