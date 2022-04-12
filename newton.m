function x = newton(f,fg,fg2_inv,x0,epsilon,epoch_max)
x = x0;

f_new = f(x);
flag = 1;
epoch = 0;
disp(['initial, f_new = ' num2str(f_new)])
while flag
    x_old = x;
    f_old = f_new;
    epoch = epoch + 1;
    
    grad = fg(x);
    step = -fg2_inv*grad;
    x = x + step;
    f_new = f(x);
    if mod(epoch,1) == 0
        disp(['epoch = ' num2str(epoch) ', f_new = ' num2str(f_new),', norm_step = ' num2str(norm(step))])
    end
    if abs(f_new-f_old)<epsilon || epoch == epoch_max
        flag = 0;
    end
end
end