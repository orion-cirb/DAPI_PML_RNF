package DAPI_PML_RNF_Tools;

import ij.ImagePlus;
import java.util.HashMap;
import mcib3d.geom2.Object3DComputation;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.image3d.ImageHandler;

/**
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class Nucleus {
    
    public Object3DInt nucleus;
    public Objects3DIntPopulation pmlFoci;
    public Objects3DIntPopulation rnf4Foci;
    
    public Nucleus(Object3DInt nucleus, Objects3DIntPopulation pmlFoci, Objects3DIntPopulation rnf4Foci) {
        this.nucleus = nucleus;
        this.pmlFoci = pmlFoci;
        this.rnf4Foci = rnf4Foci;
    }
    
    public HashMap<String, Double> computeParams(ImagePlus imgPml, ImagePlus imgRnf4, ImagePlus imgRnf111, double pixArea) {
        HashMap<String, Double> params = new HashMap<>();
        
        // Nucleus
        params.put("nucId", (double) nucleus.getLabel());
        double nucArea = new MeasureVolume(nucleus).getVolumeUnit();
        params.put("nucArea", nucArea);
        
        // PML foci
        ImageHandler imhPml = ImageHandler.wrap(imgPml);
        params.put("pmlInt", new MeasureIntensity(nucleus, imhPml).getValueMeasurement(MeasureIntensity.INTENSITY_SUM));
        params.put("pmlFociNb", (double) pmlFoci.getNbObjects());
        params.put("pmlFociArea", computePopTotArea(pmlFoci));
        params.put("pmlFociInt", computePopTotInt(pmlFoci, imhPml));
        double[] pmlDiffuseParams = computeFociDiffuseParams(nucleus, pmlFoci, imgPml);
        params.put("pmlDiffuseArea", nucArea - pmlDiffuseParams[0]);
        params.put("pmlDiffuseInt", pmlDiffuseParams[1]);
                
        // RNF4 foci
        ImageHandler imhRnf4 = ImageHandler.wrap(imgRnf4);
        params.put("rnf4Int", new MeasureIntensity(nucleus, imhRnf4).getValueMeasurement(MeasureIntensity.INTENSITY_SUM));
        params.put("rnf4FociNb", (double) rnf4Foci.getNbObjects());
        params.put("rnf4FociArea", computePopTotArea(rnf4Foci));
        params.put("rnf4FociInt", computePopTotInt(rnf4Foci, imhRnf4));
        double[] rnf4DiffuseParams = computeFociDiffuseParams(nucleus, rnf4Foci, imgRnf4);
        params.put("rnf4DiffuseArea", nucArea - rnf4DiffuseParams[0]);
        params.put("rnf4DiffuseInt", rnf4DiffuseParams[1]);
        
        // RNF111 intensity
        ImageHandler imhRnf111 = ImageHandler.wrap(imgRnf111);
        params.put("rnf111Int", new MeasureIntensity(nucleus, imhRnf111).getValueMeasurement(MeasureIntensity.INTENSITY_SUM));
        
        // PML/RNF4 foci colocalization
        colocalizePmlRnf4(pmlFoci, rnf4Foci, params, pixArea);
                   
        return(params);
    }
    
    
    /**
     * Find total volume of objects in population  
     */
    private double computePopTotArea(Objects3DIntPopulation pop) {
        double totArea = 0;
        for(Object3DInt obj: pop.getObjects3DInt()) {
            totArea += new MeasureVolume(obj).getVolumeUnit();
        }
        return(totArea);
    }
    
    
    /**
     * Find total integrated intensity of objects in population  
     */
    private double computePopTotInt(Objects3DIntPopulation pop, ImageHandler imh) {
        double totInt = 0;
        for(Object3DInt obj: pop.getObjects3DInt()) {
            totInt += new MeasureIntensity(obj, imh).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        }
        return(totInt);
    }
    
    
    /**
     * Compute raw integrated density in nucleus:
     * - Fill foci voxels with zero in foci channel
     * - Compute raw integrated density of nucleus in resulting foci channel
     */
    public double[] computeFociDiffuseParams(Object3DInt nuc, Objects3DIntPopulation fociPop, ImagePlus imgFoci) {
        ImageHandler imhFoci = ImageHandler.wrap(imgFoci.duplicate());
        double dilatedFociArea = 0;
        for (Object3DInt foci: fociPop.getObjects3DInt()) {
            Object3DInt dilatedObj = new Object3DComputation(foci).getObjectDilated(2f, 2f, 0);
            dilatedObj.drawObject(imhFoci, 0);
            dilatedFociArea += new MeasureVolume(dilatedObj).getVolumeUnit();
        }
        double fociDiffuseInt = new MeasureIntensity(nuc, imhFoci).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
        double[] params = {dilatedFociArea, fociDiffuseInt};
        imhFoci.closeImagePlus();
        return(params);
    }
    
    
        /**
     * Find partner coloc with pml
     */
    public void colocalizePmlRnf4(Objects3DIntPopulation pmlPop, Objects3DIntPopulation rnf4Pop, HashMap<String, Double> params, double pixArea) {
        double pmlNb = 0;
        double pmlArea = 0;
        double rnf4Nb = 0;
        double rnf4Area = 0;
        double colocArea = 0;
        
        MeasurePopulationColocalisation coloc = new MeasurePopulationColocalisation(pmlPop, rnf4Pop);
        for (Object3DInt pml: pmlPop.getObjects3DInt()) {
            for (Object3DInt rnf4: rnf4Pop.getObjects3DInt()) {
                double colocVal = coloc.getValueObjectsPair(pml, rnf4);
                if (colocVal > 0.25*pml.size()) {
                    pmlNb++;
                    pmlArea += new MeasureVolume(pml).getVolumeUnit();
                    break;
                }
            }
        }
        
        for (Object3DInt rnf4: rnf4Pop.getObjects3DInt()) {
            for (Object3DInt pml: pmlPop.getObjects3DInt()) {
                double colocVal = coloc.getValueObjectsPair(pml, rnf4);
                if (colocVal > 0.25*rnf4.size()) {
                    rnf4Nb++;
                    rnf4Area += new MeasureVolume(rnf4).getVolumeUnit();
                    break;
                }
            }
        }
        
        for (Object3DInt rnf4: rnf4Pop.getObjects3DInt()) {
            for (Object3DInt pml: pmlPop.getObjects3DInt()) {
                double colocVal = coloc.getValueObjectsPair(pml, rnf4);
                if ((colocVal > 0.25*pml.size()) || (colocVal > 0.25*rnf4.size())) {
                    colocArea += colocVal*pixArea;
                }
            }
        }
        
        // RNF4-positive PML foci
        params.put("pmlRnf4FociNb", pmlNb);
        params.put("pmlRnf4FociArea", pmlArea);
        // PML-positive RNF4 foci
        params.put("rnf4PmlFociNb", rnf4Nb);
        params.put("rnf4PmlFociArea", rnf4Area);
        // PML-RNF4 foci overlap
        params.put("pmlRnf4OverlapArea", colocArea);        
    }
    
}
