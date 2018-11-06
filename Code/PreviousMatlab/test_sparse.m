
signal = [1,2,3,4]
pry_medi_tans(signal)


load('signal.mat');
signal_p=zeros(1,size(signal,2));
for i=1:size(signal,2)
    if signal(i)>=0
        signal_p(i)=signal(i);
    else
        signal_p(i)=0;
    end
end
signal_m=zeros(1,size(signal,2));
for i=1:size(signal,2)
    if signal(i)<=0
        signal_m(i)=-signal(i);
    else
        signal_m(i)=0;
    end
end
alpha_p=spar_posi_decomp(signal_p);
s_p=posi_reconstr(alpha_p);
alpha_m=spar_posi_decomp(signal_m);
s_m=posi_reconstr(alpha_m);
%% the test of second generation and second generation reconstruction
%M=csvread('testing2.csv');
%wt=sec_starlet_trans(M(:,2)');
%S=posi_reconstr(wt);
%plot(10.^M(:,1)',M(:,2)');
%hold on
%plot(10.^M(:,1),S);
