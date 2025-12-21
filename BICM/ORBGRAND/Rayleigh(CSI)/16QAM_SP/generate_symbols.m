function Z = generate_symbols(angles, ams, probs, h, a, sigma, N)
    cdf = cumsum(probs);
    U = rand(1, N);
    idx = arrayfun(@(u) find(u < cdf, 1), U);
    theta = angles(idx);
    am = ams(idx);
    X = h.*am.*a.*cos(theta) + sigma * randn(1, N);
    Y = h.*am.*a.*sin(theta) + sigma * randn(1, N);
    Z = X + 1i * Y;
    Z = Z.';
end
