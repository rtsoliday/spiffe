# spiffe

**spiffe** is a 2.5-dimensional particle-in-cell (PIC) code developed at the Advanced Photon Source (APS) for the design and simulation of radio-frequency (RF) guns. This software is tailored to model the behavior of charged particles under the influence of electromagnetic fields, with a particular focus on the dynamics within RF guns. **spiffe** generates particle output files in the SDDS (Self-Describing Data Sets) format, which can be directly read and analyzed by the **[elegant](https://github.com/rtsoliday/elegant)** software for downstream beamline simulations and analysis.

## Overview
**spiffe** is designed to assist researchers and engineers in the field of accelerator physics by providing a computational tool to simulate and optimize RF gun designs. Its 2.5D approach strikes a balance between computational efficiency and the accuracy needed to capture essential physics, making it a valuable resource for studying particle dynamics in RF-driven systems.

## Installation
Detailed instructions for installing **spiffe**, including system requirements, dependencies, and compilation steps, are available in the **[spiffe manual](https://ops.aps.anl.gov/manuals/spiffe_latest/spiffe.html)**. To get started, clone the repository from GitHub:
``` bash
git clone https://github.com/rtsoliday/spiffe.git
```
Then, refer to the manual for specific guidance on building and installing the software on your system.

## Usage
**spiffe** operates by reading an input file that defines the simulation parameters, such as the RF gun geometry, electromagnetic field configurations, particle source properties, and boundary conditions. Once the simulation is executed, **spiffe** produces output files in the SDDS format, which can be post-processed using **elegant** or other SDDS-compatible tools.

A typical workflow might look like this:

1. Prepare an input file with the necessary simulation parameters.

2. Run the **spiffe** executable with the input file.

3. Analyze the resulting SDDS output files with **elegant**.

For detailed instructions on creating input files, running simulations, and interpreting outputs, please consult the **[spiffe manual](https://ops.aps.anl.gov/manuals/spiffe_latest/spiffe.html)**.

## Documentation
The  **[spiffe manual](https://ops.aps.anl.gov/manuals/spiffe_latest/spiffe.html)** is the primary resource for comprehensive information about the software. It covers:
- The theoretical background of RF gun design and PIC simulations.
- Details of the 2.5D approximation implemented in **spiffe**.
- Specifications for input file formats and parameters.
- Compilation and execution instructions.
- Guidance on analyzing SDDS output files.
Users are encouraged to review the manual for in-depth insights into **spiffe**’s capabilities and usage.

## Contributing
**spiffe** is hosted on GitHub, and contributions from the community are welcome. If you encounter bugs, have suggestions for enhancements, or wish to contribute features, please:
- Open an issue on the GitHub repository to report problems or propose ideas.
- Submit a pull request with your changes for review.

Your feedback and contributions can help improve **spiffe** for the broader accelerator physics community.

## Authors
- M. Borland
- R. Soliday
- H. Shang

## Acknowledgments
This project is developed and maintained by **[Accelerator Operations & Physics](https://www.aps.anl.gov/Accelerator-Operations-Physics)** at the **Advanced Photon Source** at **Argonne National Laboratory**.

For more details, visit the official **[SDDS documentation](https://www.aps.anl.gov/Accelerator-Operations-Physics/Documentation)**.
