syms x0 
R=300.4;
f=@(x0) (-0.5+x0)/((sqrt(R^2-x0^2))-(1-0.466)*R)-(2*x0/(sqrt(R^2-x0^2)))/(1-x0^2/(sqrt(R^2-x0^2))^2)
qvjian_start = 0;
qvjian_end = 150;
zero = [];
for chu = [qvjian_start:1:qvjian_end]
    zero = [zero,fzero(f,chu)];
end
zero = unique(zero)