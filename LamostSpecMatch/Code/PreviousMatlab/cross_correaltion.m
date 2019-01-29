% Redshift Estimation by Cross-Correlation
%input the testing data: signal S (with redshift)
%output the redshift of the signal
function z=cross_correaltion(S)% a vector
K=csvread('training2.csv');
temp=zeros(size(K,2)-1,size(K,1));
chi2=zeros(1,size(S,2));
for i=2:size(K,2)
    temp(i-1,:)=K(:,i);
end
temp=ppca(temp,20);%make a decision of N; ???
for j=1:size(S)
    for i=1:20
    b=b+(inner_prod(S,temp(:,i),j))^2;
    end
    chi2(j)=b;
end
[y,J]=max(chi2);
log_lambda=K(1,1)+J*0.0002;
z=10^(log_lambda)-1;

    
    
   









