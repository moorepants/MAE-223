function y = findzero(b, c, x0)

options = optimset('Display', 'off'); % Turn off Display
y = fzero(@poly, x0, options);

    function y = poly(x) % Compute the polynomial.
    y = x^3 + b*x + c;
    end
end
