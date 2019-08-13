##basic functions ad definitions to work with tipping problems


#definition of the TippingProblem class
struct TippingProblem
	f::Function			#f(x, lambda)
	sigma::Real			#size of the noise
	lambda::Function	#parameter shift
	lmax::Real			#limit of Lambda for t -> infty
	
	#generic constructor
	#check that sigma and lmax are > 0
	TippingProblem(F, sigma, lambda, lmax) = begin
		if sigma < 0
			error("sigma must be positive")
		elseif lmax < 0
			error("lmax must be positive")
		else
			new(f, sigma, lambda, lmax)
		end
	end
end


#drawing random samples of the attractor at time t=0
#input variables
#	F:			tipping problem
#	x0:			starting point for pullback iteration
#	npoints:	number of points to be drawn
#	nsteps:		number of pullback steps
#return variable
#	ret:		1D array of length npoints
function attractorsample(F::TippingProblem, x0::Real, npoints::Int, nsteps::Int)
	ret = fill(x0, npoints);
	g(u, x) = F.f(x, 0) + F.sigma * u;
	for n in 1:npoints, k in 1:nsteps
		ret[n] = g(2 * rand() - 1, ret[n])
	end
	return ret
end

#generate sequence of parameters Lambda(n*r) that are !=0 and != lambda_max
#input variables
#	F:			tipping problem
#	r:			rate
#return variable
#	ret:		1D array of lambda-values
function parametersequence(F::TippingProblem, r::Real)
	ret = Float64[];
	k = 0;	#running index
	
	#evlauate lambda in rn + total offset
	lbd(n) = F.lambda(r * n);
	
	#find first index with lbd(k)>0
	while lbd(k) == 0
		k += 1;
	end
	
	#until lbd(k) = lmax add lbd(k) to ret
	while lbd(k) < F.lmax
		push!(ret, lbd(k));
		k += 1;
	end
	
	#return
	return ret
end
	
#doestipalong simulates one path and decides if there is tipping or not
#input variables
#	F:			tipping problem
#	x0:			initial value for the iteration
#	par:		sequence of parameters lambda to use
#	qmin:		lower bound for repeller
#	qmax:		upper bound for repeller
#return values
#	true:		system tips
#	false:		system does not tip
function doestipalong(F::TippingProblem, x0::Real, par::Array{T,1}, qmin::Real, qmax::Real) where T <: Real
	x = x0;	#temporary variable
	
	#iterate along parameters par
	for p in par
		x = F.f(x, p) + F.sigma * (2* rand() - 1);
	end
	
	#iterate until exiting the interval
	while qmin < x < qmax
		x = F.f(x, F.lmax) + F.sigma * (2 * rand() - 1);
	end
	
	#return true iff tipping
	if x <= qmin
		return true
	else
		return false
	end
end

#calculates tipping probability for a given rate
#input variables
#	F:			tipping problem
#	attsmpl:	1D vector of sample points drawn from the attractor
#	r:			rate
#	ntrials:	number of times the experiment is run
#	qmin:		lower bound for repeller
#	qmax:		upper bound for repeller	
#return variable
#	ret:		calculated empirical tipping frequency	
function tippingprob(F, attsmpl, r, ntrials, qmin, qmax)
	
	#initialize things
	ctr = 0.0;		#counter for the number of trials performed
					#stored as float because we divide by it
	ret = 0.0;		#value to be returned, constantly updated
	#same parametersequence is used throughout
	par = parametersequence(F, r);
	#shorthand notation
	g(x) = doestipalong(F, x, par, qmin, qmax);
	
	#perform the test
	for a in attsmpl, n in 1:ntrials
		ctr += 1;	#increase the counter
		#println(a, " ", n, " ", g(a), " ", ctr, " ")
		ret += 1/ctr * (g(a) - ret)
	end
	
	#return value
	return ret
end
	








































		

