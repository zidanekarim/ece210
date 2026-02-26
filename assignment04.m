%% 1. dotprod

dotprod = @(x, y) (x)' * y;


%% 2
function [orthonormal] = is_orthonormal(array, func) 
    num_cols = size(array, 2);
    tolerance = 1000 * eps;
    for k = 1:num_cols
        if abs(sqrt(func(array(:,k), array(:,k))) - 1) > tolerance
            orthonormal = false;
            return
        end
        for j = k+1:num_cols
            if abs(func(array(:,k), array(:,j))) > tolerance 
                orthonormal = false;
                return
            end
        end
    end
    orthonormal = true;
    return
end

%% 3 gram-schmidt
function [orthonormalArray] = gram_schmidt(array, func)
    [num_rows, num_cols] = size(array);
    if (is_orthonormal(array, func))
        orthonormalArray = array;
        return
    end
    orthonormalArray = zeros(num_rows, num_cols);
    for k = 1:num_cols
        v_k = array(:, k);
        for j = 1:k-1
            v_k = v_k - (func(orthonormalArray(:, j), v_k) * orthonormalArray(:, j));
        end
        orthonormalArray(:, k) = v_k / sqrt(func(v_k, v_k)); 
    end

    
        
end

%% 4 tester

test_matrix = randi(15,4,4) + 1j * randi(15, 4 ,4);
gram_schmidt_test_matrix = gram_schmidt(test_matrix, dotprod)
output = is_orthonormal(gram_schmidt_test_matrix, dotprod)