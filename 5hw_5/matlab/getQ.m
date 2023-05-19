function Q = getQ(n_seg, n_order, ts)
    Q = [];
    for k = 1:n_seg
        Q_k = zeros(n_order+1,n_order+1);
        %#####################################################
        % STEP 1.1: calculate Q_k of the k-th segment 
        %
        %
        %
        %
        for i=4:n_order
            for l=4:n_order
                I=factorial(i)/factorial(i-4);
                L=factorial(l)/factorial(l-4);
                q_il=I*L*(ts(k)^(i+l-7))/(i+l-7);
                Q_k(i+1,l+1)=q_il;
            end
        end
        Q = blkdiag(Q, Q_k);
    end
end