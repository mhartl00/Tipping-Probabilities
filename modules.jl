#VARYING SIGMA
#module for calculating graphs with fixed lambda_max and several values of sigma
module VarySigma
#get definitions and parameters
include("functions.jl")
include("parameters.jl")


#define the sigmas
sigma_min = 0;
sigma_step = 0.01;
sigma_length = 6;
if sigma_min + sigma_step * sigma_length > sigma_max
	print("some sigmas exceed sigma_max")
end
sigma = range(sigma_min, step = sigma_step, length = sigma_length);


#calculate tipping probability as function of rate
#input variables
#	lmax:	fixed value of lambda_max to use
#	test:	use test parameters (test = true) or regular parameters (test = false)
#output variable: 2D array
#	first line: lmax, sigma[1], ..., sigma[n]
#	from second line on 
#		- first columns has rates
#		- all other columns hold the calculated probabilities
function varysigma(lmax; test = false)
	#set parameters
	plot_parameters(test);
	plot_rates(test);

	#generate tipping problems
	F = Array{Any, 1}(undef, sigma_length);
	for n in 1:sigma_length
		F[n] = TippingProblem(f, sigma[n], get_lambda(lmax), lmax);
	end

	#set up the return array
	p = Array{Float64, 2}(undef, nrates + 1, sigma_length + 1)
	p[1, 1] = lmax
	p[1, 2:end] = sigma
	p[2:end, 1] = rates
	
	#do the calculations
	for n in 1:sigma_length
		smpl = attractorsample(F[n], 1.0, npoints, nsteps);	#here the sample depends on sigma
		prob(r) = tippingprob(F[n], smpl, r, ntrials, qmin(F[n]), qmax(F[n]))
		p[2:end, n + 1] = prob.(rates);
	end

	#return value
	return p
	
end #end of varysigma()

end #end of module VarySigma





#VARYING LAMBDA_MAX
#module for calculating graphs with fixed sigma and several values of lambda_max
module VaryLmax

#get definitions and parameters
include("functions.jl")
include("parameters.jl")

#parameters
sigma = 0.04;
if sigma > sigma_max
	print("sigma is too large")
end

lmax_min = 1;
lmax_step = 0.05;
lmax_length = 6;


#construct tipping problems
lmax = range(lmax_min, step = lmax_step, length = lmax_length);
F = Array{Any, 1}(undef, lmax_length);
for n in 1:lmax_length
	F[n] = TippingProblem(f, sigma, get_lambda(lmax[n]), lmax[n]);
end


#calculate tipping probability as function of rate
#input variables
#	sigma:	fixed value of sigma to use
#	test:	use test parameters (test = true) or regular parameters (test = false)
#output variable: 2D array
#	first line: sigma, lambda[1], ..., lambda[n]
#	from second line on 
#		- first columns has rates
#		- all other columns hold the calculated probabilities
function varylmax(; test = false)
	#set plot parameters
	plot_parameters(test)
	plot_rates(test)
	
	#attractor sample
	#NB: does not depend on lmax
	spl = attractorsample(F[1], 1.0, npoints, nsteps)
	
	#set up the return array
	p = Array{Float64, 2}(undef, nrates + 1, lmax_length + 1)
	p[1, 1] = sigma
	p[1, 2:end] = lmax
	p[2:end, 1] = rates
	
	#do the calculations
	for n in 1:lmax_length
		#shortcut for tipping probability function
		prob(r) = tippingprob(F[n], spl, r, ntrials, qmin(F[n]), qmax(F[n]))
		#apply to rates
		p[2:end, n + 1] = prob.(rates)
	end
	
	#return value
	return p
end	#end of varylmax()

end #end of module VaryLmax



	













