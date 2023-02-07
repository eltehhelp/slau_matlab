function [x, ok] = lab_slau_kramer(A, b, tol)
    [m, n] = size(A);
    ok = true;
    x = zeros(m, 1);
    
    if m ~= n || abs(det(A)) < tol
        ok = false;
       
            fprintf('ERROR: Kramers method is not applicable.\n\n');
        disp(x)
    else
        x = zeros(size(b));

        for i = 1 : 1 : m
            M = A;
            M(:, i) = b;
            x(i) = det(M) / det(A);
        end

        %if show
         %   fprintf('Решение для метода Крамера:\n');
          %  disp(x);
        %end
    end


end

