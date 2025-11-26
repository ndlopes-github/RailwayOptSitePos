# Release Notes

## Version 0.5.0 - November 26, 2025

[![DOI](https://zenodo.org/badge/659697901.svg)](https://doi.org/10.5281/zenodo.17725507)

This release accompanies the publication of the research paper in the 4OR journal.

### Paper Publication

**Title:** Minimizing costs in signal provision from communication antennas along a railway line

**DOI:** [10.1007/s10288-025-00599-7](https://doi.org/10.1007/s10288-025-00599-7)

**Authors:** A. Araújo, J. O. Cerdeira, N. Lopes, A. Moura

**Journal:** 4OR - A Quarterly Journal of Operations Research

**Publication Date:** September 24, 2025

### Overview

This release provides a complete implementation of the algorithm described in the paper, including:

- 0/1 linear optimization model for wireless network design on railway lines
- Real data examples (Section 3.1)
- Simulated data generation and experiments (Section 3.2)
- Visualization tools for signal coverage analysis

### Key Features

#### Real Data Implementation
- Model execution for real-world railway line instances
- Signal coverage visualization with RX signal level plots
- Fair and good coverage interval projections
- CSV output for optimal facility selection
- Incremental summary tables

#### Simulated Data Generation
- Seed-based reproducible data generation
- Configurable instance characteristics
- Batch simulation capabilities
- Automated parameter control

#### Output and Reporting
- Console logging and summaries
- Signal coverage graphs
- Gurobi solver logs
- CSV data files for analysis
- Summary tables for paper reproduction

### Technical Requirements

- Julia Language
- DrWatson package for project management
- Gurobi optimizer (or HiGHS as open-source alternative)
- Bash shell for batch scripts

### Installation

```bash
git clone https://github.com/ndlopes-github/RailwayOptSitePos.git
cd RailwayOptSitePos
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Usage

**Real Data Examples:**
```bash
julia scripts/model_solvit.jl
```

**Simulated Data Generation:**
```bash
bash scripts/simsdatagen.sh
bash scripts/loopsims.sh
```

### Known Limitations

- Data generation and simulations are computationally intensive
- Default scripts use reduced instance counts for testing purposes
- Active Gurobi license required (HiGHS alternative available)

### File Structure

```
RailwayOptSitePos/
├── scripts/          # Execution scripts
├── data/             # Input and output data
│   ├── exp_pro/      # Generated instances
│   └── sims/         # Simulation results
├── plots/            # Signal coverage visualizations
└── src/              # Core implementation
```

### Documentation

Usage instructions and examples are available in the [README.md](README.md).

### Citation

If you use this code in your research, please cite:

```bibtex
@article{araujo2025minimizing,
  title={Minimizing costs in signal provision from communication antennas along a railway line},
  author={Ara{\'u}jo, A. and Cerdeira, J. O. and Lopes, N. and Moura, A.},
  journal={4OR},
  year={2025},
  doi={10.1007/s10288-025-00599-7}
}
```

### Contributing

Contributions are welcome! Please open an issue or submit a pull request.

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Contact

For questions or inquiries: nuno(dot)lopes(at)isel(dot)pt

---

**Note:** This is a stable release accompanying the published research. Future updates may include performance improvements and additional features.