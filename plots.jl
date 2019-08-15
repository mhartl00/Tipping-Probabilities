#routines for plotting the graphs
using Plots

using DelimitedFiles
using LaTeXStrings

include("modules.jl")
import .VaryLmax

pyplot()

#setting the default xticks
default_xtick_pos = [n/6 for n in 0:4]
default_xtick_labels = [
	L"0",
	L"\frac{1}{6}",
	L"\frac{1}{3}",
	L"\frac{1}{2}",
	L"\frac{2}{3}"
]
default_xticks = (default_xtick_pos, default_xtick_labels)

#basic plotting function
function pr_plot(; 
	rates::AbstractArray{T}, 
	pr::AbstractArray{T}, 
	title::AbstractString, 
	label,
	xlabel::AbstractString = L"r", 
	ylabel::AbstractString = L"\mathrm{Pr}\,(r)",
	xticks = default_xticks,
) where T <: Real
	return plot(
		rates, 
		pr,
		windowsize = (800, 400),
		#aspect_ratio = .5,
		label = label,
		title = title,
		xlabel = xlabel,
		ylabel = ylabel,
		xticks = xticks,
		xtickfont = font(13),
		ytickfont = font(11),
		xlim = [rates[1], rates[end]],
		ylim = [0, 1],
		margin = 7Plots.mm,
		reuse = false,
		show = true
	)
end


#make the plot for varying lambda max
using .VaryLmax

function vary_lmax_plot(p)
	#preprocess data
	lmax = p[1, 2:end]
	n = length(lmax)
	sigma = p[1, 1]
	rates = p[2:end, 1]
	pr = p[2:end, 2:end]
	
	#make plotting labels
	labs = Array{Any, 2}(undef, 1, n)
	for n in 1:n
		labs[n] = latexstring("\\lambda_+ =", lmax[n])
	end
	
	#make the plot
	plt = pr_plot(
		rates = rates,
		pr = pr,
		label = labs,
		title = latexstring("\\sigma =", sigma)
	)
end


#plot for varying sigma
using .VarySigma

function vary_sigma_plot(p)
	#preprocess data
	sigma = p[1, 2:end]
	n = length(sigma)
	lmax = p[1, 1]
	rates = p[2:end, 1]
	pr = p[2:end, 2:end]
	
	#make plotting labels
	labs = Array{Any, 2}(undef, 1, n)
	for n in 1:n
		labs[n] = latexstring("\\sigma =", sigma[n])
	end
	
	#make the plot
	plt = pr_plot(
		rates = rates,
		pr = pr,
		label = labs,
		title = latexstring("\\lambda_+ =", lmax)
	)
end
