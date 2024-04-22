import DAPI_PML_RNF_Tools.Nucleus;
import DAPI_PML_RNF_Tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.DebugTools;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/**
 * Detect DAPI nuclei
 * Detect PML and RNF4 foci
 * Compute colocalization between PML and RNF4 foci
 * @author Héloïse Monnet @ ORION-CIRB
 */
public class DAPI_PML_RNF implements PlugIn {
    
    DAPI_PML_RNF_Tools.Tools tools = new Tools();

    public void run(String arg) {
        try {  
            if (!tools.checkInstalledModules() || !tools.checkStardistModels()) {
                return;
            }
                        
            String imageDir = IJ.getDirectory("Choose images directory")+File.separator;
            if (imageDir == null) {
                IJ.showMessage("Specified directory not found");
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(new File(imageDir));
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            DebugTools.setRootLevel("warn");
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata(); 
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] chNames = tools.dialog(channels);
            if (chNames == null) {
                IJ.showMessage("Plugin canceled");
                return;
            }
            
            // Create output folder for results files and images
            String outDir = imageDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write headers in results files
            BufferedWriter results = new BufferedWriter(new FileWriter(outDir + "results.csv", false));
            results.write("Image name\tNucleus ID\tNucleus area (µm2)\tPML raw integrated density\tPML foci number\tPML foci total area (µm2)\t"
                    + "PML foci raw integrated density\tPML diffuse area (µm2)\tPML diffuse raw integrated density\t"
                    + "RNF4 raw integrated density\tRNF4 foci number\tRNF4 foci total area (µm2)\tRNF4 foci raw integrated density\t"
                    + "RNF4 diffuse area (µm2)\tRNF4 diffuse raw integrated density\tRNF111 raw integrated density\t"
                    + "RNF4-positive PML foci number\tRNF4-positive PML foci total area (µm2)\tPML-positive RNF4 foci number\t"
                    + "PML-positive RNF4 foci total area (µm2)\tPML-RNF4 foci total overlap area (µm2)\n");
            results.flush();
            
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                System.out.println("--- ANALYZING IMAGE " + rootName + " ---");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Open all channels
                tools.print("- Opening RNF4 channel -");
                int index = ArrayUtils.indexOf(channels, chNames[0]);
                ImagePlus imgRnf4 = BF.openImagePlus(options)[index];
                tools.print("- Opening PML channel -");
                index = ArrayUtils.indexOf(channels, chNames[1]);
                ImagePlus imgPml = BF.openImagePlus(options)[index];
                tools.print("- Opening RNF111 channel -");
                index = ArrayUtils.indexOf(channels, chNames[2]);
                ImagePlus imgRnf111 = BF.openImagePlus(options)[index];
                tools.print("- Opening DAPI channel -");
                index = ArrayUtils.indexOf(channels, chNames[3]);
                ImagePlus imgDapi = BF.openImagePlus(options)[index];
                
                // Detect DAPI nuclei
                tools.print("- Detecting DAPI nuclei -");
                Objects3DIntPopulation nucPop = tools.cellposeDetection(imgDapi);
                
                // Detect PML foci
                tools.print("- Detecting PML foci -");
                Objects3DIntPopulation pmlPop = tools.fociDetection(imgPml);
                
                // Detect RNF4 foci
                tools.print("- Detecting RNF4 foci -");
                Objects3DIntPopulation rnf4Pop = tools.stardistDetection(imgRnf4);
                
                // Get nuclei with their respective population of PML and RNF4 foci
                tools.print("- Getting PML and RNF4 foci for each nucleus -");
                ArrayList<Nucleus> nuclei = tools.colocalizeNucFoci(nucPop, pmlPop, rnf4Pop);
               
                // Write results
                tools.print("- Writing results -");
                for (Nucleus nucleus: nuclei) {
                    HashMap<String, Double> p = nucleus.computeParams(imgPml, imgRnf4, imgRnf111, tools.pixArea);
                    results.write(rootName+"\t"+p.get("nucId").intValue()+"\t"+p.get("nucArea")+"\t"+p.get("pmlInt")+"\t"+p.get("pmlFociNb").intValue()+"\t"+p.get("pmlFociArea")+
                            "\t"+p.get("pmlFociInt")+"\t"+p.get("pmlDiffuseArea")+"\t"+p.get("pmlDiffuseInt")+"\t"+p.get("rnf4Int")+"\t"+p.get("rnf4FociNb").intValue()+"\t"+p.get("rnf4FociArea")+
                            "\t"+p.get("rnf4FociInt")+"\t"+p.get("rnf4DiffuseArea")+"\t"+p.get("rnf4DiffuseInt")+"\t"+p.get("rnf111Int")+"\t"+p.get("pmlRnf4FociNb").intValue()+
                            "\t"+p.get("pmlRnf4FociArea")+"\t"+p.get("rnf4PmlFociNb").intValue()+"\t"+p.get("rnf4PmlFociArea")+"\t"+p.get("pmlRnf4OverlapArea")+"\n");
                    results.flush();
                }
                
                // Draw results
                tools.print("- Drawing results -");
                tools.drawResults(nuclei, imgRnf4, imgPml, imgDapi, imgRnf111, outDir, rootName);
                
                tools.closeImage(imgRnf4);
                tools.closeImage(imgPml);
                tools.closeImage(imgRnf111);
                tools.closeImage(imgDapi);
            }
            System.out.println("--- All done! ---");
        } catch (NullPointerException | IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(DAPI_PML_RNF.class.getName()).log(Level.SEVERE, null, ex);
            tools.print("Plugin aborted");
        }
    }
}
