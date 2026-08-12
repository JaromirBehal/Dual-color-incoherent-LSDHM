```markdown
# Data and Function Structure

The repository is organized as follows:

## Main MATLAB file

**`MultiWavelengths_Cheek_cells_and_Oil_drop.m`**

The main MATLAB script used for loading and processing the camera images. It performs the reconstruction and analysis of the experimental data.

## `data/`

Contains the raw camera images acquired during the experiments with cheek cells and the oil drop.

## `functions/`

Contains the MATLAB functions used by the main script:

- **`apodiz_SG.m`** — function for apodization of the input data.
- **`apodization_for_propag.m`** — function for applying apodization prior to numerical propagation.
- **`off_axis_reconstruction_general.m`** — general function for off-axis digital holographic reconstruction.
