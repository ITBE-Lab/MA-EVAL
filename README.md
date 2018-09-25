# <img src="https://raw.githubusercontent.com/ITBE-Lab/MA/release/MA.png" align="center" width="90"> The Modular Aligner Evaluation Tool

This is a collection of code snippets that can be used to analyze aligners. This project needs to 
be configured in the *settings.py* file. Further, the paths in *command_line_aligner.py* need to be 
adjusted (it is possible to add or remove aligners by adding/deleting further child classes for
CommandLine).

## Brief summary of file contents
- *create_accuracy_graph.py* generates the accuracy and runtime analysis.
- *misc.py* computes the startup times and memory usage.
- *ambiguity_distrib.py* computes seed ambiguity for a genome. (requires MA with WITH_PYTHON=1)
- *SV_analysis.py* evaluates SV mapping quality.

## Underlying data for figures in *Accurate high throughput alignment via line sweep-based seed processing*
The underlying data for all figures in **link to final publication** can be downloaded here: _link_

