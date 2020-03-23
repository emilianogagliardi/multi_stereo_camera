function [val,jacob] = objective_function(x)
    global FUN;
    global JACOBIAN;
    val = FUN(x);
    if nargout > 1
        jacob = JACOBIAN(x);
    end
end

