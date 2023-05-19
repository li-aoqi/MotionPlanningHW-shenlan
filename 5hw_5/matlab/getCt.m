function Ct = getCt(n_seg, n_order)
    %#####################################################
    % STEP 2.1: finish the expression of Ct
    %
    r=(n_order+1)*n_seg;
    c= 8+(n_seg-1)*4;
    Ct = zeros(r,c);
    
    for i=1:4
        Ct(i,i)=1;
    end
    
    for i=1:4
        Ct(8*(n_seg-1)+4+i,4+(n_seg-1)+i)=1;
    end
    
    if n_seg~=1
        for i=0:n_seg-2

            Ct(8*i+4+1,4+i+1)=1;
            Ct(8*(i+1)+1,4+i+1)=1;

            Ct(8*i+4+2,8+(n_seg-1)+i*3+1)=1;
            Ct(8*(i+1)+2,8+(n_seg-1)+i*3+1)=1;

            Ct(8*i+4+3,8+(n_seg-1)+i*3+2)=1;
            Ct(8*(i+1)+3,8+(n_seg-1)+i*3+2)=1;

            Ct(8*i+4+4,8+(n_seg-1)+i*3+3)=1;
            Ct(8*(i+1)+4,8+(n_seg-1)+i*3+3)=1;
        end
    end
end