% 求解贝塞尔多项式系数
function poly_coef = MinimumSnapCorridorBezierSolver(axis, waypoints, corridor, ts, n_seg, n_order, v_max, a_max)
    % 起点和终点条件
    start_cond = [waypoints(1), 0, 0];
    end_cond   = [waypoints(end), 0, 0];   
    
    %% #####################################################
    % STEP 1: 计算对应于普通多项式的Q矩阵和对应于贝塞尔曲线的Q_0矩阵
    [Q, M]  = getQM(n_seg, n_order, ts);
    Q_0 = M'*Q*M;
    Q_0 = nearestSPD(Q_0);
    
    %% #####################################################
    % STEP 2: 获取相应的等式约束
    [Aeq, beq] = getAbeq(n_seg, n_order, ts, start_cond, end_cond);
    
    %% #####################################################
    % STEP 3: 获取飞行走廊中各个区域的范围以及相应的不等式约束
    
    % STEP 3.1: get corridor_range of x-axis or y-axis,
    % you can define corridor_range as [p1_min, p1_max;
    %                                   p2_min, p2_max;
    %                                   ...,
    %                                   pn_min, pn_max ];
    corridor_range = [];
    for k = 1:n_seg
        corridor_range(k,:) = [corridor(axis,k) - corridor(axis+2,k), corridor(axis,k) + corridor(axis+2,k)];
    end
    % STEP 3.2: 获取相应的不等式约束
    [Aieq, bieq] = getAbieq(n_seg, n_order, corridor_range, ts, v_max, a_max);
    
    f = zeros(size(Q_0,1),1);
    poly_coef = quadprog(Q_0,f,Aieq, bieq, Aeq, beq);
end