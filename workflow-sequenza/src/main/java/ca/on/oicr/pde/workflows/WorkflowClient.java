package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
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
    private String outputFileNamePrefix;

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
    private String bin;
    private String pypy;
    private String rScript;
    private String rLib;

    //ref Data
    private String refFasta;
    private String sequenzaGCData;

    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(WorkflowClient.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;

    // meta-types
    private final static String PDF_METATYPE = "application/pdf";
    private final static String TXT_METATYPE = "text/plain";
    private final static String TAR_GZ_METATYPE = "application/tar-gzip";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            tumorBam = getProperty("input_files_tumor");
            normalBam = getProperty("input_files_normal");

            //Ext id
            outputFileNamePrefix = getProperty("output_file_prefix");

            //bin data 
            sequenzaGCData = getProperty("sequenza_bin_data_hg19");

            //samtools
            samtools = getProperty("samtools");

            // ref fasta
            refFasta = getProperty("ref_fasta");

            //r path
            rScript = getProperty("rpath") + "/bin/Rscript";
            rLib = getProperty("rLib");
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
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("tumor");
        file0.setSourcePath(tumorBam);
        file0.setType("application/bam");
        file0.setIsInput(true);
        SqwFile file1 = this.createFile("normal");
        file1.setSourcePath(normalBam);
        file1.setType("application/bam");
        file1.setIsInput(true);
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        /**
         * Steps for sequenza:
         * 1. Check if "bam" file exists; true
         * 2. Check if "bai" file exists; true: go to step 4
         * 3. Check if normal Pb_R sample exists; true: go to step 4; else abort
         * 3. If false: samtools index "bam" file
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
        String inputTumorBamFilePath = getFiles().get("tumor").getProvisionedPath();
        String inputNormalBamFilePath = getFiles().get("normal").getProvisionedPath();
        String externalIdentifier = this.outputFileNamePrefix;
        String outputDir = externalIdentifier + "_output";
        String tempDir = this.tmpDir;

        String sample_name = inputTumorBamFilePath.substring(inputTumorBamFilePath.lastIndexOf("/") + 1, inputTumorBamFilePath.lastIndexOf(".bam"));
        intermediateFilePath = tempDir + sample_name + "seqz.bin50.gz";

        Job sequenzaUtilJob = getSequenzaUtilsJob(intermediateFilePath, inputTumorBamFilePath, inputNormalBamFilePath);
        //sequenzaUtilJob.addParent(parentJob);
        parentJob = sequenzaUtilJob;

        Job runSequenzaR = runSequenzaRJob(intermediateFilePath, outputDir);
        runSequenzaR.addParent(parentJob);
        parentJob = runSequenzaR;

        Job zipOutput = iterOutputDir(outputDir);
        zipOutput.addParent(parentJob);

        // Provision *.pdf, .seg, model-fit.tar.gz files
        String segFile = this.outputFileNamePrefix + "_Total_CN.seg";
        String pdfFile = this.outputFileNamePrefix + "_genome_view.pdf";
        SqwFile cnSegFile = createOutputFile(outputDir + "/" + segFile, TXT_METATYPE, this.manualOutput);
        cnSegFile.getAnnotations().put("segment data from the tool ", "Sequenza ");
        zipOutput.addFile(cnSegFile);

        SqwFile cnImage = createOutputFile(outputDir + "/" + pdfFile, PDF_METATYPE, this.manualOutput);
        cnImage.getAnnotations().put("copy number view ", "Sequenza ");
        zipOutput.addFile(cnImage);

        SqwFile zipFile = createOutputFile(outputDir + "/" + "model-fit.tar.gz", TAR_GZ_METATYPE, this.manualOutput);
        zipFile.getAnnotations().put("Other files ", "Sequenza ");
        zipOutput.addFile(zipFile);
    }

    // create Job function for the sequenza steps - pre-step
    private Job getSequenzaUtilsJob(String intFilePath, String inputTumorBamFilePath, String inputNormalBamFilePath) {
        Job jobSequenzaUtils = getWorkflow().createBashJob("sequenza-utils");
        Command command = jobSequenzaUtils.getCommand();
        command.addArgument(pypy);
        command.addArgument(sequenzaUtil);
        command.addArgument("bam2seqz");
        command.addArgument("--fasta " + refFasta);
        command.addArgument("-n " + inputNormalBamFilePath);
        command.addArgument("-t " + inputTumorBamFilePath);
        command.addArgument("-gc " + sequenzaGCData);
        command.addArgument("-S " + samtools + " |");
        command.addArgument(pypy);
        command.addArgument(sequenzaUtil);
        command.addArgument("seqz-binning");
        command.addArgument("-s -");
        command.addArgument("-w 50 | gzip > " + intFilePath);
        jobSequenzaUtils.setMaxMemory(Integer.toString(sequenzaUtilMem * 1024));
        jobSequenzaUtils.setQueue(getOptionalProperty("queue", ""));
        return jobSequenzaUtils;
    }

    private Job runSequenzaRJob(String intFilePath, String outDir) {
        Job jobSequenzaR = getWorkflow().createBashJob("sequenza_R");
        Command cmd = jobSequenzaR.getCommand();
        cmd.addArgument("export R_LIBS=" + rLib + ";");
        cmd.addArgument(rScript);
        cmd.addArgument(sequenzaRscript);
        cmd.addArgument(outDir);
        cmd.addArgument(intFilePath);
        cmd.addArgument(outputFileNamePrefix);
        jobSequenzaR.setMaxMemory(Integer.toString(sequenzaRscriptMem * 1024));
        jobSequenzaR.setQueue(getOptionalProperty("queue", ""));
        return jobSequenzaR;
    }

    private Job iterOutputDir(String outDir) {
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
        Job iterOutput = getWorkflow().createBashJob("handle_output");
        Command cmd = iterOutput.getCommand();
        cmd.addArgument("bash -x " + getWorkflowBaseDir() + "/dependencies/handleFile.sh");
        cmd.addArgument(this.outputFileNamePrefix);
        cmd.addArgument(outDir);
        iterOutput.setMaxMemory(Integer.toString(sequenzaRscriptMem * 1024));
        iterOutput.setQueue(getOptionalProperty("queue", ""));
        return iterOutput;
    }

}
