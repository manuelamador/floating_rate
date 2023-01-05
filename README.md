# Sovereign debt crises and floating-rate bonds 

The repository contains the code associated with the paper:

["Sovereign debt crises and floating-rate bonds"](https://manuelamador.me/files/floatingrate.pdf)
    
by Mark Aguiar, Manuel Amador and Ricardo Alves Monteiro (2023).


## Structure

The subfolder `src` contains the main source code.

The subfolder `scripts` contains some of the analysis of the model for certain parameters. It contains both a julia script as well as the corresponding jupyter notebook. The scripts generate the figures and moments reported in the paper.  

The subfolder `output` contains the figures generates by the scripts, as well as some the calculated moments. 

## Running the code 

The code is in [Julia](https://julialang.org/downloads/).

To run the code, open a julia prompt at the root of this repository and type:

    julia> using Pkg 
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()

The above will download the packages needed to run the code. 
  

To run the jupyter notebook, do:
  
    julia> using IJulia
    julia> notebook(dir=".")
  
That should open a browser with [Jupyter](https://jupyter.org/) . Navigate to `scripts` to locate the notebooks. 

There are five notebooks in `scripts`:

  - `floating_rate_figures_tables.ipynb` solves the one period bond model (with and without runs), the long term bond model (with and without runs) and the floating rate bond model. It generates the plots and tables presented in  ["Sovereign debt crises and floating-rate bonds"](https://manuelamador.me/files/floatingrate.pdf)

   

