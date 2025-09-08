# RailwayOptSitePos

This repository reflects the implementation of an algorithm based on a paper that is set to appear. Please note that the contents of this repository are currently under development and may undergo changes as we refine the implementation to align with the paper. We encourage you to keep an eye on this repository for updates.

## Paper Details

   > Title: Minimizing costs in signal provision from communication antennas along a railway line

   

   > Authors: A. Araújo 1, J. O. Cerdeira 2, N. Lopes 3, A. Moura 4

   > 1- CMUC, Department of Mathematics, University of Coimbra;<br>
   > 2- CMA, Department of Mathematics, NOVA University Lisbon;<br>
   > 3- ISEL, Polytechnic of Lisbon, and CEMAT, University of Lisbon;<br>
   > 4- ISEP-LEMA, Polytechnic of Porto, and CMUP, University of Porto;

   > Journal: [4OR](https://link.springer.com/journal/10288)

   > Publication Date: [To appear]

   > Abstract:
In this paper we address a wireless network design problem on a railway line. Given a finite set of locations along a railway line and different types of communication antennas that can be installed at each of these locations, which locations and which type of antenna should be selected to ensure a certain level of signal coverage along the railway line while minimizing construction costs?
We formulate the problem as a 0/1 linear optimization model, prove that the problem is NP-hard, and report computational experiments using real and simulated data. The computational tests showed that the model is capable of solving the problem for railway lines longer than any existing real railway lines. 

# Description

The code in this repository aims to replicate the algorithm described in the forthcoming paper. Our intention is to provide a practical implementation that can be used and tested by the community. As such, please consider this code as a work in progress, subject to further modifications and improvements.


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

- Note: An active Gurobi license is required to run the code as is. However, you can replace Gurobi with HiGHS, an open-source optimization solver. To do this, install the HiGHS.jl package and modify the solver settings in the scripts accordingly. Refer to the HiGHS documentation for detailed instructions.

# USAGE

##  Real Data (see Section 3.1 of the paper):

### Running the Model

To run the model, follow the steps below:

0. Open a console or terminal at the project directory.

1. Run the following command:
   ```
   $julia scripts/model_solvit.jl
   ```
      >   This command executes instance 4 of the model with $LMAX^n=16.2$ and $LMIN^g=0$ (see Table 2 of the paper for more details).

* Note: The previous command can be replaced by the usage of an IDE such as VSCode or by running the script in the Julia REPL.

### Output

After running the model, you can expect the following output:

+ A summary of the results and log information displayed at the terminal or console.

+ Graphs representing the signal coverage of all facilities and the  signal coverage provided by the optimal solutions saved in the plots directory.

   Example: Signal coverage of the facility selection

   ![RX Signal Level of the antennas of the facility selection](plots/solution_real_data_signal.png)

     ![Intervals of fair and good coverage of the facility selection](plots/solution_real_data_projection.png)




+ A CSV file reflecting the optimal facility selection  saved at  data/sims/solvit directory.

+ An incremental table  with the summary of the optimal solutions saved at data/sims/solvit/table.txt and data/sims/solvit/raw_table.csv.

## Simulated Data (see Section 3.2 of the paper):

### Generating and Simulating Data

To reproduce the examples in Section 3.2, follow these steps:

0. Open a console or terminal at the project directory.

1. **Generate Data**

   Run the following script to generate the simulated data:
   ```
   $ bash scripts/simsdatagen.sh
   ```
   > This script calls the Julia script `sim_instances_generator.jl`.  
   > The randomization process is seed-based to ensure reproducibility.

   The generated data will be located in:
   ```
   data/exp_pro/t1_12_8-t2_15_10-p_4_2_nD
   ```

2. **Generate Simulations**

   Run the following script to generate the simulations:
   ```
   $ bash scripts/loopsims.sh
   ```
   > This script calls the Julia script `model_sim.sh`.

### Output

After running the simulations, you can expect the following output:

+ A summary of the results and log information displayed at the terminal or console.

+ Tables and Gurobi log files saved in:
  ```
  data/sims/t1_12_8-t2_15_10-p_4_2_nD
  ```
 (These outputs are used to generate the correspondent tables in the paper.)

**Additional Notes**

- The number of generated instances is controlled by the `simsdatagen.sh` script.
- The characteristics of the generated data are controlled by the `sim_instances_generator.jl` script.
- The number of simulations is controlled by the `loopsims.sh` script.
- The model parameters are controlled by the `model_sim.sh` script.
- To facilitate testing, the number of generated instances in the scripts has been reduced, as data generation and simulations are time- and resource-intensive. Adjust these values in the scripts as needed for your experiments.

# Contributing

We welcome contributions to this project! If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request. We appreciate your involvement in making this implementation more robust and accurate.

# License

This project is currently under MIT license. Please refer to the LICENSE file for more information.

# Contact

If you have any questions or inquiries regarding this codebase or the associated paper, please contact: [nuno(dot)lopes(at)isel(dot)pt]
