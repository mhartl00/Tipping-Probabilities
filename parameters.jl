using Roots
using Optim

#setting the plot parameters
function plot_parameters(test = false)
	if test
		global npoints = 100;
		global nsteps = 100;
		global ntrials = 100;
	else 		
		global npoints = 2000; 
		global nsteps = 1000; 
		global ntrials = 5000;
	end
end

#setting the rates along which to evaluate tipping probabilities
function plot_rates(test = false)
	if test
		global nrates = 100;
	else
		global nrates = 450;
	end
	global rates = range(0.01, stop = 0.7, length = nrates);
end

#one step function f of our example
arccot(x) = pi/2 - atan(x);
ptilde(x) = x * arccot(x) + .5 * log(1 + x^2) + .1 * x;
pp(x) = ptilde(x) / ptilde(1)
f(x, p) = pp(x - p) + p

#parameter shift Lambda as a function of lambda_max
function get_lambda(lmax)
	return function(t)
		if t < 1/3
			return 0
		elseif t > 2/3
			return lmax
		else
			return lmax * (3 * t - 1)
		end
	end
end

#calculate sigma_max for the example
res = optimize(x -> x - f(x,0), 0, 1, GoldenSection());	#finds the minimum of -(f(x, 0) - x)
sigma_max = - Optim.minimum(res)						#getting the result and make it the maximum

#calculate q_min and q_max for example
qmin(F::TippingProblem) = find_zero(x -> F.f(x, 0) - x + F.sigma, 0) + F.lmax
qmax(F::TippingProblem) = find_zero(x -> F.f(x, 0) - x - F.sigma, 0) + F.lmax
