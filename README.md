# RailwayOptSitePos

This repository reflects the implementation of an algorithm based on a paper that is set to appear. Please note that the contents of this repository are currently under development and may undergo changes as we refine the implementation to align with the paper. We encourage you to keep an eye on this repository for updates.

## Paper Details

   > Title: Minimizing costs in signal provision from communication antennas along a railway line

   > Authors: A. Araújo 1, J. O. Cerdeira 2, N. Lopes 3, A. Moura 4
   > 1- CMUC, Department of Mathematics, University of Coimbra;
    2- CMA, Department of Mathematics, NOVA University Lisbon;
    3- ISEL, Polytechnic of Lisboa, and CEMAT, University of Lisboa;
    4- ISEP-LEMA, Polytechnic of Porto, and CMUP, University of Porto;
   
   > Journal: [To appear]
   > Publication Date: [To appear]

# Description

The code in this repository aims to replicate the algorithm described in the forthcoming paper. Our intention is to provide a practical implementation that can be utilized and tested by the community. As such, please consider this code as a work in progress, subject to further modifications and improvements.


# Installation

## Julia and DrWatson
This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> RailwayOptSitePos

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "RailwayOptSitePos"
```
which auto-activate the project and enable local path handling from DrWatson.

# USAGE 

## Real World Data (see Section 3.1):
### TO APPEAR
      julia scripts/model_solvit.jl

##  Simulated Data (see Section 3.2):

### Running the Model

To run the model, follow the steps below:

0. Open a console or terminal at the project directory.
1. Create the directories: mkdir plots data data/sims data/exp_pro

      > These directories are required in order to save the output files.

2. Run the following command: julia scripts/model_sim.jl

   >  This command executes one instance of the model with n = 64 (refer to Table 3 for more details).

* Note: The previous command can be replaced by the usage of an IDE such as VSCode or by running the script in the Julia REPL.

### Output

After running the model, you can expect the following output:

+ A summary of the results will be displayed in the terminal or console.

+ Graphs representing the signal coverage of all facilities and the the signal coverage provided by the optimal solutions, are be saved in the plots directory.
   
   Example: Signal coverage of the facility selection
   ![Signal coverage of the facility selection](aux/solution_coverage.png)

+ A CSV file reflecting the selected solution will be saved in the data/sims directory.

+ The file data/exp_pro/table.txt will be incremented with a summary of the solution and processing times.


# Contributing

We welcome contributions to this project! If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request. We appreciate your involvement in making this implementation more robust and accurate.
License

This project is currently under [choose license] license. Please refer to the LICENSE file for more information.

# Contact

If you have any questions or inquiries regarding this codebase or the associated paper, please contact [Insert contact information].
