function div_res = alpha_beta_div(x,y,alpha,beta)
    % Check that the parameters of the alpha-beta divergence are real numbers
    if ~isreal(alpha) || ~isreal(beta)
        my_stack = dbstack;
        fct_name = my_stack.name;
        error("Error in %s: The parameters alpha and beta must be real scalar numbers.", fct_name);
    end
      
    % Handle different cases of alpha and beta
    if alpha == 0 && beta == 0
        div_res = 0.5 * sum((log(x) - log(y)).^2, 'all'); % Log-Euclidean distance
    elseif alpha == 0 && beta ~= 0
        div_res = 1 / beta^2 * sum(y.^beta .* log(y.^beta ./ x.^beta) - y.^beta + x.^beta, 'all'); % Beta divergence with alpha = 0
    elseif alpha ~= 0 && beta == 0
        div_res = 1 / alpha^2 * sum(x.^alpha .* log(x.^alpha ./ y.^alpha) - x.^alpha + y.^alpha, 'all'); % Beta divergence with beta = 0
    elseif alpha == -beta && alpha ~= 0
        div_res = 1 / alpha^2 * sum(log(y.^alpha ./ x.^alpha) + (y.^alpha ./ x.^alpha) - 1, 'all'); % Special case for alpha = -beta
    elseif alpha ~= 0 && beta ~= 0 && alpha + beta ~= 0
        div_res = -1 / (alpha * beta) * sum(x.^alpha .* y.^beta - alpha / (alpha + beta) * x.^(alpha + beta) - beta / (alpha + beta) * y.^(alpha + beta), 'all'); % General case
    end
end


% function div_res = alpha_beta_div(x, y, alpha, beta)
%     % 计算 beta 散度
%     if alpha == 1 && beta == 1
%         div_res = sum(x .* log(x ./ y) - x + y, 'all'); % KL 散度
%     elseif alpha == 0 && beta == 0
%         div_res = 0.5 * sum((x - y).^2, 'all'); % 欧几里得距离的平方
%     elseif alpha == 0 && beta == 1
%         div_res = sum(x - y + y .* log(y ./ x), 'all'); % Itakura-Saito 散度
%     elseif alpha == 1 && beta == 0
%         div_res = sum(x ./ y - log(x ./ y) - 1, 'all'); % Itakura-Saito 散度
%     else
%         div_res = -1 / (alpha * beta) * sum(x.^alpha .* y.^beta - alpha / (alpha + beta) * x.^(alpha + beta) - beta / (alpha + beta) * y.^(alpha + beta), 'all'); % General case
%     end
% end
