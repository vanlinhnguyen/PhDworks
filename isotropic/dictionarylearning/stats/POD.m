function [D, A, lambda] = POD(X)
    C = (X'*X)./(size(X,1)-1);                  %'# cov(X)

    [D, lambda] = eig(C);
    [lambda, order] = sort(diag(lambda), 'descend');       %# sort cols high to low
    D = D(:,order);

    A = X*D(:,1:end);
end