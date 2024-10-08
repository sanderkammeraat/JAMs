Simple proof-of-concept simulations of Active Matter in Julia.  

**Disclaimer:** This code is a hobby project and one of the goals is to learn the Julia programming language, so it likely is not fully optimized (yet).  

See examples on how to use.  

**Possible workflow**  
See https://code.visualstudio.com/docs/languages/julia for detailed instructions  
1) Install Julia.  
2) Install Visual Studio Code.  
3) Install Julia extension for VS Code.
4) Install Julia packages (e.g. via REPL terminal in VS Code).  
5) For proper live plotting, disable the "Use Plot Pane" option in the settings of the Julia extension of VS code.


**Required Julia packages:**  
-StaticArrays (for speedup for small vectors, e.g. position or force vectors)  
-Plots (for live plotting)  
-ProgressBars (for progressbar in integrator)  
-Random (for noise and random initial positions, probably pre-installed)  
-Distributions (for noise and random initial positions, probably pre-installed)  
-LinearAlgebra (for doing linear algebra, probably pre-installed)  

**Multi-threading**  
For information about multi-threading in Julia and how to set the number of threads, visit:  
https://docs.julialang.org/en/v1/manual/multi-threading/  


