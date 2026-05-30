# JAMs: Simple proof-of-concept simulations of Active Matter in Julia.

<img width="1398" height="1376" alt="sim_snapshot" src="https://github.com/user-attachments/assets/11490a11-a9f0-419b-9787-5ffa8049ed38" />

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
-GLMakie (for plotting)  
-Observables (for live plotting)  
-Distributions (for noise and random initial positions)  
-JLD2 to save simulation settings for loading simulations within Julia  and to use JLD2 features when loading the raw simulation data.   
-HDF5 (for outputting raw simulation data and metadata in HDF5 format).  
-CodecZlib (for saving with compression)  
-ProgressMeter (for providing a progressbar, ETA and it/s)  

The packages mentioned above can be installed in one go by using the REPL command:  
`add StaticArrays GLMakie Observables Distributions JLD2 HDF5 CodecZlib ProgressMeter`  

Probably pre-installed:  
-LinearAlgebra (for doing linear algebra)  
-Random (for noise and random initial positions)  


**Multi-threading**  
For information about multi-threading in Julia and how to set the number of threads, visit:  
https://docs.julialang.org/en/v1/manual/multi-threading/  


