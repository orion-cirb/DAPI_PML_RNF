# DAPI_PML_RNF

* **Developed for:** Hsin-Chieh
* **Team:** De Thé
* **Date:** April 2024
* **Software:** Fiji

### Images description

2D images of primary cells taken with a x60 objective

4 channels:
  1. *Alexa Fluor 488:* RNF4 foci
  2. *Alexa Fluor 594:* PML foci
  3. *Alexa Fluor 647:* RNF111
  4. *DAPI:* DAPI nuclei

### Plugin description

* Detect DAPI nuclei using Cellpose
* Detect PML foci using Difference of Gaussians filtering + Otsu thresholding + Fill Holes filtering
* Detect RNF4 dots using Stardist
* Compute colocalization between PML and RNF4 foci
* For each nucleus, provide area and intensity in various channels, foci area and intensity in corresponding channel, diffuse area and intensity, and overlap area between PML and RNF4 foci 

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ** Fiji plugin
* **Stardist** model named *pmls2.zip*
* **Cellpose** conda environment + *cyto* pretrained model
   

### Version history

Version 1 released on April 15, 2024.
