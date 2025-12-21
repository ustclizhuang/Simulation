function Z = generate_symbols(angles, probs, h, a, sigma, N)
    cdf = cumsum(probs);
    U = rand(1, N);
    idx = arrayfun(@(u) find(u < cdf, 1), U);
    theta = angles(idx);
    X = sqrt(2)*h*a*cos(theta) + sigma * randn(1, N);
    Y = sqrt(2)*h*a*sin(theta) + sigma * randn(1, N);
    Z = X + 1i * Y;
    Z = Z.'; 
end
