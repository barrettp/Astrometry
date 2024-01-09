all:
	julia --project=. -e 'import Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
	julia --project=. make.jl

clean:
	rm -rf build/

preview:
	julia -e 'using LiveServer; serve(dir="build")'
