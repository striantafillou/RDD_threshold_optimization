function yint =y_int(x, a, b, h)
    % estimate utility for linear function 
    tmp1 = b*h^2+2*a*h;
    tmp2 = b.*(x.^2)+2*a.*x;
    yint = 0.5*(tmp1-tmp2);
end