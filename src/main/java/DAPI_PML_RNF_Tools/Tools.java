package DAPI_PML_RNF_Tools;

import DAPI_PML_RNF_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import DAPI_PML_RNF_Tools.Cellpose.CellposeTaskSettings;
import DAPI_PML_RNF_Tools.StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class Tools {
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/DAPI_PML_RNF";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    String[] chNames = {"RNF4: ", "PML: ",  "RNF111: ", "DAPI: "};
    public Calibration cal;
    public double pixArea;
    
    // Nuclei detection with Cellpose
    private String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    private final String cellposeModelPath = IJ.isWindows()? System.getProperty("user.home")+"\\.cellpose\\models\\" : "";
    private String cellposeModel = "cyto";
    private int cellposeDiam = 350;
    private double nucAreaMin = 30;
    private double nucAreaMax = 150;
    private double nucIntMin = 5000;
    
    // PML foci detection with DoG + Thresholding
    private final double dogSigma1A = 1;
    private final double dogSigma1B = 3;
    private final double dogSigma2A = 3;
    private final double dogSigma2B = 6;
    private final String thMethod = "Otsu";
    
    // RNF4 foci detection with Stardist
    private final double scalingFactor = 0.75;
    private final Object syncObject = new Object();
    private final File stardistModelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    private final String stardistOutput = "Label Image"; 
    private final String stardistModel = "pmls2.zip";
    private final double stardistPercentileBottom = 0.2;
    private final double stardistPercentileTop = 99.8;
    private double stardistProbThresh = 0.8;
    private double stardistOverlapThresh = 0.25;
    
    // Foci size filtering
    private double fociAreaMin = 0.01;
    private double fociAreaMax = 1;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom2.Object3DInt");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Check that required StarDist models are present in Fiji models folder
     */
    public boolean checkStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = stardistModelsPath.listFiles(filter);
        int index = ArrayUtils.indexOf(modelList, new File(stardistModelsPath+File.separator+stardistModel));
        if (index == -1) {
            IJ.showMessage("Error", stardistModel + " StarDist model not found, please add it in Fiji models folder");
            return false;
        }
        return true;
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
 
    
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth);
        return(cal);
    }
    
    
    /**
     * Find channels name
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 60, 0);
        gd.addImage(icon);

        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[n]);
        }
        
        gd.addMessage("Nuclei detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min nucleus area (µm2): ", nucAreaMin, 2);
        gd.addNumericField("Max nucleus area (µm2): ", nucAreaMax, 2);
        gd.addNumericField("Min nucleus intensity: ", nucIntMin, 2);
        
        gd.addMessage("Foci detection", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min foci area (µm2): ", fociAreaMin);
        gd.addNumericField("Max foci area (µm2): ", fociAreaMax);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm): ", cal.pixelHeight, 4);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chChoices = new String[chNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();

        nucAreaMin = gd.getNextNumber();
        nucAreaMax = gd.getNextNumber();
        nucIntMin = gd.getNextNumber();
        
        fociAreaMin = gd.getNextNumber();
        fociAreaMax = gd.getNextNumber();

        cal.pixelHeight = cal.pixelWidth = gd.getNextNumber();
        pixArea = cal.pixelHeight*cal.pixelWidth;
        
        if (gd.wasCanceled())
            return(null);
        return(chChoices);
    }
    
    
    /**
     * Detect objects in 2D using Cellpose
     */
    public Objects3DIntPopulation cellposeDetection(ImagePlus imgIn) {
        // Define CellPose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModelPath+cellposeModel, 1, cellposeDiam, cellposeEnvDir);
        settings.useGpu(true);

        // Run Cellpose
        ImagePlus img = imgIn.duplicate();
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
        ImagePlus imgOut = cellpose.run();        
        imgOut.setCalibration(cal);
        
        // Filter objects 
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        System.out.println("Nb nuclei detected: "+pop.getNbObjects());
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(nucAreaMin/pixArea, nucAreaMax/pixArea);
        pop = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(img), false);
        System.out.println("Nb nuclei remaining after filtering: "+ pop.getNbObjects());
        pop.resetLabels();
        
        closeImage(img);
        closeImage(imgOut);
        return(pop);
    }
    
    
    /**
     * Detect dots with DoG + thresholding
     */
    public Objects3DIntPopulation fociDetection(ImagePlus imgIn) {
        ImagePlus imgDOG1 = DOG(imgIn, dogSigma1A, dogSigma1B);
        ImagePlus imgDOG2 = DOG(imgIn, dogSigma2A, dogSigma2B);
        ImagePlus imgDOG = ImageCalculator.run(imgDOG1, imgDOG2, "Max create");
        ImagePlus imgBin = threshold(imgDOG, thMethod);
        ImagePlus imgFill = fillHoles(imgBin);
        imgFill.setCalibration(cal);
        
        Objects3DIntPopulation pop = getPopFromImage(imgFill);
        System.out.println("Nb PML foci detected: "+pop.getNbObjects());
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(fociAreaMin/pixArea, fociAreaMax/pixArea);
        System.out.println("Nb PML foci remaining after filtering: "+ pop.getNbObjects());
        pop.resetLabels();
        
        closeImage(imgDOG1);
        closeImage(imgDOG2);
        closeImage(imgDOG);
        closeImage(imgBin);
        closeImage(imgFill);
        return(pop);
    }

    
    /**
     * Difference of Gaussians filtering using CLIJ2
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCL);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    /**
     * Fill holes using CLIJ2
     */
    public ImagePlus fillHoles(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.binaryFillHoles(imgCL, imgCLBin);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
        
    /**
     * Return population of 3D objects population from binary image
     */
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        ImageInt labels = new ImageLabeller().getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        labels.closeImagePlus();
        return(pop);
    }
    
    
    /**
     * Detect objects in 2D using Stardist
     */
    public Objects3DIntPopulation stardistDetection(ImagePlus img) throws NullPointerException {
        // Resize image to be in a StarDist-friendly scale
        ImagePlus imgIn = img.resize((int)(img.getWidth()*scalingFactor), (int)(img.getHeight()*scalingFactor), 1, "average");
        IJ.run(imgIn, "Add Slice", "add=slice");// TODO: correct bug with StardistOrion when single image is provided
           
        // Define Stardist settings
        File starDistModelFile = new File(stardistModelsPath+File.separator+stardistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgIn);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThresh, stardistOverlapThresh, stardistOutput);
        
        // Run Stardist
        star.run();
        ImagePlus imgLabels = star.getLabelImagePlus().resize(img.getWidth(), img.getHeight(), 1, "none");
        imgLabels.setSlice(2);
        IJ.run(imgLabels, "Delete Slice", "");
        imgLabels = medianFilter(imgLabels, 1);
        imgLabels.setCalibration(cal);
        
        // Filter objects 
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));
        System.out.println("Nb RNF4 foci detected: "+pop.getNbObjects());
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(fociAreaMin/pixArea, fociAreaMax/pixArea);
        System.out.println("Nb RNF4 foci remaining after filtering: "+ pop.getNbObjects());
        pop.resetLabels();
        
        closeImage(imgIn);
        closeImage(imgLabels);   
        return(pop);
    }
    
    
    /**
     * 2D median filter using CLIJ2
     */ 
    public ImagePlus medianFilter(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median2DBox(imgCL, imgCLMed, sizeXY, sizeXY);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCL);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    
    /**
     * Return list of Nucleus with their respective population of PML and RNF4 foci
     */
    public ArrayList<Nucleus> colocalizeNucFoci(Objects3DIntPopulation nucPop, Objects3DIntPopulation pmlPop, Objects3DIntPopulation rnf4Pop) {
        ArrayList<Nucleus> nuclei = new ArrayList();
        nucPop.getObjects3DInt().stream().forEach(nuc -> {
            Objects3DIntPopulation pmlInNucPop = new Objects3DIntPopulation();
            pmlPop.getObjects3DInt().stream()
                .filter(pml -> nuc.contains(new MeasureCentroid(pml).getCentroidRoundedAsVoxelInt()))
                .forEach(pml -> {pmlInNucPop.addObject(pml);});
                
            Objects3DIntPopulation rnf4InNucPop = new Objects3DIntPopulation();
            rnf4Pop.getObjects3DInt().stream()
                .filter(rnf4 -> nuc.contains(new MeasureCentroid(rnf4).getCentroidRoundedAsVoxelInt()))
                .forEach(rnf4 -> {rnf4InNucPop.addObject(rnf4);});
                
            nuclei.add(new Nucleus(nuc, pmlInNucPop, rnf4InNucPop));
        });
        return(nuclei);
    }
    
    
    /**
     * Draw results
     */
    public void drawResults(ArrayList<Nucleus> nuclei, ImagePlus imgRnf4, ImagePlus imgPml, ImagePlus imgDapi, ImagePlus imgRnf111, String outDir, String imgName) {
        ImageHandler imhNuc = ImageHandler.wrap(imgPml).createSameDimensions();
        ImageHandler imhPml = imhNuc.createSameDimensions();
        ImageHandler imhRnf4 = imhNuc.createSameDimensions();
        
        for (Nucleus nucleus: nuclei) {
            nucleus.nucleus.drawObject(imhNuc);
            for (Object3DInt pml: nucleus.pmlFoci.getObjects3DInt())
                pml.drawObject(imhPml, 255);
            for (Object3DInt rnf4: nucleus.rnf4Foci.getObjects3DInt())
                rnf4.drawObject(imhRnf4, 255);
        }
        
        imhNuc.getImagePlus().setDisplayRange(0, 1);
        ImagePlus[] imgColors = {imhRnf4.getImagePlus(), imhPml.getImagePlus(), imhNuc.getImagePlus(), imgRnf4, imgPml, imgDapi, imgRnf111};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        new FileSaver(imgObjects).saveAsTiff(outDir + imgName + ".tif");
        
        imhNuc.closeImagePlus();
        imhPml.closeImagePlus();
        imhRnf4.closeImagePlus();
        closeImage(imgObjects);
    }
        
}
    
