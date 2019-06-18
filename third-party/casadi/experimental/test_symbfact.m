ncol = 5;
nrow = 5;
colind = [0, 3, 6, 8, 10, 12];
row = [0, 1, 4, 1, 2, 4, 0, 2, 0, 3, 3, 4];
nz = [19, 12, 12, 21, 12, 12, 21, 16, 21, 5, 21, 18];
S = casadi.Sparsity(nrow, ncol, colind, row);
A = casadi.DM(S, nz);
for ata=[false, true]
    % Make s.p.d.?
    if ata
        A_test = A;
    else
        A_test = A'*A;
    end 
    S_test = A_test.sparsity();
    
    % MATLAB builtin
    if ata
        [count,h,parent,post, L] = symbfact(sparse(A_test), 'col', 'lower');
    else
        [count,h,parent,post, L] = symbfact(sparse(A_test), 'sym', 'lower');
    end
    % CasADi
    [ca_count, ca_parent, ca_post, ca_L] = S_test.symbfact(ata);
    disp('parent')
    disp(parent)
    disp(ca_parent+1)

    disp('post')
    disp(post)
    disp(ca_post+1)

    disp('count')
    disp(count)
    disp(ca_count)

    if ~ca_L.is_null() 
      disp(full(casadi.DM(ca_L, 1)))
      disp(full(L))
    end
end
