function coeff = func_pear(X,Y)  
    fenzi = sum(X.* Y) - (sum(X) * sum(Y)) / length(X);  
    fenmu = sqrt((sum(X.^2) - sum(X)^2 / length(X)) * (sum(Y.^2) - sum(Y)^2 / length(X)));  
    coeff = fenzi / fenmu; 
end