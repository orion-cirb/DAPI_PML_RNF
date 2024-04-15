# Vglut2_GluR2_GFP_Fiber

* **Developed for:** Maëla
* **Team:** Selimi
* **Date:** September 2023
* **Software:** Fiji

### Images description

3D images taken with a x60 objective

3 channels:
  1. *CSU_488:* VGlut2 dots
  2. *CSU_568:* GluR2 dots
  2. *CSU_647:* GFP neuronal fibers

### Plugin description

* Detect GFP fibers and VGlut2 dots with Quantile Based Normalization + Weka
* Detect GluR2 dots with Stardist
* Colocalize VGlut2 dots with GFP fibers, keep only those being GFP+
* Associate each VGlut2 dot with all GluR2 dots closer than 200nm (one GluR2 can only be associated with one VGlut2)

### Dependencies

* **3DImageSuite** Fiji plugin
* **Stardist** model named *fociRNA-1.2.model*
* 2 **Weka pixel classifier** models named *classifier-GFP.model* and *classifier-VGlut2.model*
   

### Version history

Version 1 released on September 5, 2023.
