function M = getM(n_seg, n_order, ts)
    M = [];
    for k = 1:n_seg
        M_k = [];
        %#####################################################
        % STEP 1.1: calculate M_k of the k-th segment 
        %
        for i=1:4
            for j=i:n_order+1
                M_k(i,j)=(factorial(j-1)/factorial(j-i))*0^(j-i);
                M_k(i+4,j)=(factorial(j-1)/factorial(j-i))*ts(k)^(j-i);
            end
        end
        M = blkdiag(M, M_k);
    end
end