function [Aeq beq]= getAbeq(n_seg, n_order, waypoints, ts, start_cond, end_cond)
    n_all_poly = n_seg*(n_order+1);
    %#####################################################
    % p,v,a,j constraint in start, 
    Aeq_start = zeros(4, n_all_poly);
    beq_start = zeros(4, 1);
    % STEP 2.1: write expression of Aeq_start and beq_start
    
    for k=0:3
        Aeq_start(k+1,k+1)=factorial(k);
    end
    beq_start=start_cond';
    %#####################################################
    % p,v,a constraint in end
    Aeq_end = zeros(4, n_all_poly);
    beq_end = zeros(4, 1);
    % STEP 2.2: write expression of Aeq_end and beq_end

    for i = 0:3
        for l=i:n_order
            Aeq_end(i+1,(n_seg-1)*(n_order+1)+l+1)= (factorial(l)/ factorial(l-i))*(ts(n_seg)^(l-i));
        end
    end
    beq_end=end_cond';

    %#####################################################
    % position constrain in all middle waypoints
    Aeq_wp = zeros(n_seg-1, n_all_poly);
    beq_wp = zeros(n_seg-1, 1);
    % STEP 2.3: write expression of Aeq_wp and beq_wp
    %
    for i=0:n_seg-2
        for l=0:n_order
            Aeq_wp(i+1,i*(n_order+1)+l+1)=(ts(i+1)^l);
        end
        beq_wp(i+1)=waypoints(i+2);
    end
    
    %#####################################################
    % position continuity constrain between each 2 segments
    Aeq_con_p = zeros(n_seg-1, n_all_poly);
    beq_con_p = zeros(n_seg-1, 1);
    % STEP 2.4: write expression of Aeq_con_p and beq_con_p
    %
    for i=0:n_seg-2
        for l=0:n_order
            Aeq_con_p(i+1,i*(n_order+1)+l+1)= ts(i+1)^l;
            Aeq_con_p(i+1,(i+1)*(n_order+1)+l+1) = -0^l;
        end
    end
    %#####################################################
    % velocity continuity constrain between each 2 segments
    Aeq_con_v = zeros(n_seg-1, n_all_poly);
    beq_con_v = zeros(n_seg-1, 1);
    % STEP 2.5: write expression of Aeq_con_v and beq_con_v
    %
    for i=0:n_seg-2
        for l=1:n_order
            Aeq_con_v(i+1,i*(n_order+1)+l+1)= l*ts(i+1)^(l-1);
            Aeq_con_v(i+1,(i+1)*(n_order+1)+l+1) = -l*0^(l-1);
        end
    end
    %#####################################################
    % acceleration continuity constrain between each 2 segments
    Aeq_con_a = zeros(n_seg-1, n_all_poly);
    beq_con_a = zeros(n_seg-1, 1);
    % STEP 2.6: write expression of Aeq_con_a and beq_con_a
    %
    for i=0:n_seg-2
        for l=2:n_order
            Aeq_con_a(i+1,i*(n_order+1)+l+1)= l*(l-1)*ts(i+1)^(l-2);
            Aeq_con_a(i+1,(i+1)*(n_order+1)+l+1) = -l*(l-1)*0^(l-2);
        end
    end
    %#####################################################
    % jerk continuity constrain between each 2 segments
    Aeq_con_j = zeros(n_seg-1, n_all_poly);
    beq_con_j = zeros(n_seg-1, 1);
    % STEP 2.7: write expression of Aeq_con_j and beq_con_j
    %
    for i=0:n_seg-2
        for l=3:n_order
            Aeq_con_j(i+1,i*(n_order+1)+l+1)= l*(l-1)*(l-2)*ts(i+1)^(l-3);
            Aeq_con_j(i+1,(i+1)*(n_order+1)+l+1) = -l*(l-1)*(l-2)*0^(l-3);
        end
    end

    %#####################################################
    % combine all components to form Aeq and beq   
    Aeq_con = [Aeq_con_p; Aeq_con_v; Aeq_con_a; Aeq_con_j];
    beq_con = [beq_con_p; beq_con_v; beq_con_a; beq_con_j];
    Aeq = [Aeq_start; Aeq_end; Aeq_wp; Aeq_con];
    beq = [beq_start; beq_end; beq_wp; beq_con];