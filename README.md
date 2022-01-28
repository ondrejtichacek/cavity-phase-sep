
### Installation

```
conda create -n cavity-phase-sep
conda activate cavity-phase-sep
conda install -c conda-forge python numpy scipy cython matplotlib coloredlogs verboselogs

pip install -e .
```

### Example usage
```
python optimize_cy.py --plot --molecule Lys10_ADP --exp-data-file sep_intensity_Poly-Lys10_ADP_4.txt
```
