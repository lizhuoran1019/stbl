function [alpha_estimate,beta_estimate, gamma_estimate, delta_estimate] = stblfmpe(impulseNoise, p)
%分数阶矩估计参数
Ex_p = mean(abs(impulseNoise).^p, 2);
q = -p;
Ex_q = mean(abs(impulseNoise).^q, 2);
left = 2*tan(p*pi/2) ./ (p*pi*Ex_p.*Ex_q);
syms a;
alpha_estimate = zeros(size(impulseNoise,1), 1);
for i = 1:size(impulseNoise,1)
    alpha_estimate(i) = double(abs( vpasolve(sin(p*pi/a) / (p*pi/a) == left(i), a) ));
end
Cpa = 2.^(p+1).*gamma((p+1)./2).*gamma(-p./alpha_estimate) ./ ( alpha_estimate.*sqrt(pi).*gamma(-p/2) );
gamma_estimate = double((Ex_p ./ Cpa).^(1./p));
beta_estimate = 0;
delta_estimate = 0;
end