% find the position where the start observation point is the same 
% just by comparing very lambda_tr with the first entry of the lambda_te
function k=get_k(lambda_tr,lambda_te)
k=1;
for k=1:size(lambda_tr,2)
    if lambda_tr(k)==lambda_te(1)
        fprintf('%d',k);
        break;
    else
        k=k+1;
    end
end
