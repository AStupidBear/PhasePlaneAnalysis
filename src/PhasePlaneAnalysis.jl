module PhasePlaneAnalysis
	using MATLAB
	cd(joinpath(Pkg.dir("PhasePlaneAnalysis"),"src"))
	spawn(`jupyter-notebook PhasePlaneAnalysis.ipynb`)
	mat"addpath($(pwd()))"
end # module
