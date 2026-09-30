# Spheroid Colonization

MATLAB scripts for quantifying bacterial colonization of cancer spheroids from fluorescence microscopy Z-stack images.

## Analysis Workflow

The analysis is performed as a four-step pipeline:

### 1. `Segment_spheroids.m`

Segments spheroid boundaries from the nuclear fluorescence channel.

**Input:**

* Multi-channel Z-stack TIFF image

**Output:**

* Spheroid segmentation masks and boundaries (`.mat`)

---

### 2. `Analyze_noninfected_controls.m`

Quantifies GFP and mCherry fluorescence in non-infected spheroids to determine background fluorescence for each Z-slice.

**Input:**

* Non-infected control TIFF image
* Spheroid segmentation mask

**Output:**

* Per-Z-slice GFP and mCherry background intensities (`.mat`)

---

### 3. `Quantify_average_BG.m`

Calculates the average background fluorescence profile from multiple non-infected control spheroids.

**Input:**

* Quantification results from non-infected controls

**Output:**

* Average GFP and mCherry background intensity for each Z-slice

---

### 4. `Quantify_colonization.m`

Quantifies bacterial fluorescence in infected spheroids and subtracts the background fluorescence measured from non-infected controls.

The script calculates:

* Background-corrected GFP fluorescence
* Background-corrected mCherry fluorescence
* Total integrated GFP and mCherry fluorescence
* GFP/mCherry fluorescence ratio

**Input:**

* Infected spheroid TIFF image
* Spheroid segmentation mask
* Averaged background fluorescence

**Output:**

* Background-corrected colonization measurements
* Summary results

## Requirements

* MATLAB
* Fluorescence microscopy Z-stack images in TIFF format
* Nuclear fluorescence channel for spheroid segmentation
* GFP and mCherry fluorescence channels for bacterial quantification

## Notes

Run the scripts in the order shown above. The output of each step is used as an input for the subsequent analysis step.

## License

This scripts are licensed under CC BY-NC-SA 4.0
