package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.io.File;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
/**
 * <p>For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare Java Workflows</a>.</p>
 * 
 * Quick reference for the order of methods called:
 * 1. setupDirectory
 * 2. setupFiles
 * 3. setupWorkflow
 * 4. setupEnvironment
 * 5. buildWorkflow
 * 
 * See the SeqWare API for 
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a> 
 * for more information.
 */
public class WorkflowClient extends OicrWorkflow {
    //dir
    private String dataDir, tmpDir;
    
    // Input Data
    private String tumorBam;
    private String normalBam;
    private String externalId;
    
    // Output check
    private boolean isFolder = true;
    
    //Scripts 
    private String sequenzaUtil;
    private String sequenzaRscript;
    
    //Tools
    private String samtools;
    
    //Memory allocation
    private Integer sequenzaUtilMem;
    private Integer sequenzaRscriptMem;
    
    //path to bin
    private String pypy;
    private String rScript;
    
    //ref Data
    private String refFasta;
    private String sequenzaGCData; 
    
    
    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(WorkflowClient.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;
    
            

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");
            
            // input samples 
            tumorBam = getProperty("input_files_tumor");
            normalBam = getProperty("input_files_tumor");

            //Ext id
            externalId = getProperty("external_name");

            //bin data 
            sequenzaGCData = getProperty("sequenza_bin_data_hg19");

            //samtools
            samtools = getProperty("samtools");

            // ref fasta
            refFasta = getProperty("ref_fasta");

            //r path
            rScript = getProperty("rpath") + "/bin/Rscript";
            pypy = getProperty("pypy");
            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");
            
            // sequenza
            sequenzaUtil = getProperty("sequenza_utils_script");
            sequenzaRscript = getProperty("sequenza_rscript");
            sequenzaUtilMem = Integer.parseInt(getProperty("sequenza_utils_mem"));
            sequenzaRscriptMem = Integer.parseInt(getProperty("sequenza_rscript_mem"));
            
            

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }
    
   
    @Override
    public void buildWorkflow() {
        
        /**
         * Steps for sequenza:
         * 1. Check if "bam" file exists; true 
         * 2. Check if "bai" file exists; true: go to step 4
         * 3. Check if normal Pb_R sample exists; true: go to step 4; else abort
         * 3. If false:  samtools index "bam" file
         * 4. Run job sequenza-utils
         * 5. If outputFile ends with "bin50.gz"; go to step 6; else go to step 4
         * 6. Run job sequenzaR
         * 7. Iterate through the files/folders in outDir:
         * 8. If fileName1 == "pandc.txt" and fileName2 ends with "Total_CN.seg"; create a folder called "copynumber"
         * 9. If fileType == "folder"; create a folder called "model-fit"; move folders to "model-fit"
         * 10. If fileType == "file" && fileName != outputFile; move file to "model-fit"
         * 11. Delete outputFile (rm outputFile)
         * 12. zip "model-fit"
         * 13. outputFile = fileName2
         * 14. OutputDir contains the following: fileName1, outputFile, model-fit.zip
         */
        
        
        // workflow : read inputs tumor and normal bam; run sequenza-utils; write the output to temp directory; 
        // run sequenzaR; handle output; provision files (3) -- model-fit.zip; text/plain; text/plain
        
        Job parentJob = null;
        String intermediateFilePath;
        String inputNormalBamFilePath = this.normalBam;
        String inputTumorBamFilePath = this.tumorBam;
        String externalIdentifier = this.externalId;
        String outputDir = externalIdentifier+"_output";
        String tempDir = this.tmpDir;
        
        String sample_name = inputTumorBamFilePath.substring(inputTumorBamFilePath.lastIndexOf("/") + 1, inputTumorBamFilePath.lastIndexOf(".bam"));
        intermediateFilePath = tempDir+sample_name + "seqz.bin50.gz";
        
        Job sequenzaUtilJob = getSequenzaUtilsJob(inputNormalBamFilePath, inputTumorBamFilePath, intermediateFilePath);
        sequenzaUtilJob.addParent(parentJob);
        parentJob = sequenzaUtilJob;
        
        Job runSequenzaR = runSequenzaRJob(intermediateFilePath, outputDir);
        runSequenzaR.addParent(parentJob);
        parentJob = runSequenzaR;
        
        String cmd = iterOutputDir(outputDir);
        Job zipFiles = getWorkflow().createBashJob("zip-model-fit");
        zipFiles.addParent(parentJob);
        Command command = zipFiles.getCommand();
        command.addArgument(cmd);

    }
    
    // create Job function for the sequenza steps - pre-step
    private Job getSequenzaUtilsJob(String normalSampleBamFilePath, String tumorSampleBamFilePath, String intFilePath) {
        Job jobSequenzaUtils = getWorkflow().createBashJob("sequenza-utils");
        Command command = jobSequenzaUtils.getCommand();
        command.addArgument(pypy);
        command.addArgument(sequenzaUtil);
        command.addArgument("bam2seqz");
        command.addArgument("--fasta " + refFasta);
        command.addArgument("-n " + normalSampleBamFilePath);
        command.addArgument("-t " + tumorSampleBamFilePath);
        command.addArgument("-gc " + sequenzaGCData);
        command.addArgument("-S " + samtools + " |");
        command.addArgument(pypy);
        command.addArgument(sequenzaUtil);
        command.addArgument("seqz-binning");
        command.addArgument("-s");
        command.addArgument("- -w 50 | gzip > " + intFilePath);
        jobSequenzaUtils.setMaxMemory(Integer.toString(sequenzaUtilMem*1024));
        jobSequenzaUtils.setQueue(getOptionalProperty("queue", ""));
        return jobSequenzaUtils;
    }

    private Job runSequenzaRJob(String intFilePath, String outDir) {
        Job jobSequenzaR = getWorkflow().createBashJob("sequenza_R");
        Command cmd = jobSequenzaR.getCommand();
        cmd.addArgument(rScript);
        cmd.addArgument(sequenzaRscript);
        cmd.addArgument(outDir);
        cmd.addArgument(intFilePath);
        cmd.addArgument(externalId);
        jobSequenzaR.setMaxMemory(Integer.toString(sequenzaRscriptMem*1024));
        jobSequenzaR.setQueue(getOptionalProperty("queue", ""));
        return jobSequenzaR;
    }
    
    String iterOutputDir(String outDir) {
        /**
         * Method to handle file from the output directory All provision files
         * are in tempDir Create a directory called model-fit in output
         * directory 
         * move the subfolders into it 
         * move files with the following
         * extentions to model-fit
         * "_log.txt", ".pdf", "_solutions.txt", "_CP.txt",
         * "_mutations.txt", "segments.txt", ".RData"
         * construct a cmd string to zip the model-fit folder
         */
        // find only folders in the output Directory
        String[] items = {"_log.txt", ".pdf", "_solutions.txt", "_CP.txt", "_mutations.txt", "segments.txt", ".RData"};

        File modelDir = new File(outDir + "/model-fit/");
        boolean foo = modelDir.mkdir();
        File dir = new File(outDir);
        File[] directoryListing = dir.listFiles();
        if (directoryListing != null) {
            for (File sub : directoryListing) {
                if (sub.isDirectory()) {
                    modelDir = new File(modelDir.toString() + sub.toString());
                    if (foo) {
                        sub.renameTo(modelDir);
                    }
                } else {
                    for (String item : items) {
                        if (sub.toString().endsWith(item)) {
                            modelDir = new File(modelDir.toString() + sub.toString());
                            if (foo) {
                                sub.renameTo(modelDir);
                            }
                        }
                    }
                }
            }
        }
        String modelFitDir = outDir + "/model-fit/";
        //File zipFolder = new File(modelFitDir);
        //File[] listFolder = zipFolder.listFiles();
        String cmd = null;
        //if (listFolder != null) {
        cmd = "zip -r ";
        cmd += modelFitDir;
        return cmd;
    }

}
