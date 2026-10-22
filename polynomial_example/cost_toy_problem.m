function y = cost_toy_problem(a,data,x_domain)
    model_output = a(1) + a(2)*x_domain + a(3)*x_domain.^2 + a(4)*x_domain.^7;
    residual = (model_output - data).^2;
    y = 0.5*sum(residual);
end