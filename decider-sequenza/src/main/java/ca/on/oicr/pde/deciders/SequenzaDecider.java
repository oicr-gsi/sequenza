package ca.on.oicr.pde.deciders;

import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;
import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;

/**
 *
 * @author prath@oicr.on.ca
 */
public class SequenzaDecider extends OicrDecider {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private String templateType = "EX";
    private String queue = "";
    private String sampleName;
    private Boolean groupByAligner = true;

    private final static String BAM_METATYPE = "application/bam";
    private String tumorType;
    private List<String> duplicates;
    private String tumourFilePath;
    private String normalFilePath;
    private String commaSeparatedFilePaths;
    private String commaSeparatedParentAccessions;
    private String intervalBed;
    private String refGenome = "/.mounts/labs/PDE/data/gatkAnnotationResources/hg19_random.fa";
    private String rsConfigXMLPath = "/.mounts/labs/PDE/data/rsconfig.xml";
    private Rsconfig rs;

    public SequenzaDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("template-type", "Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type. Default: " + this.templateType).withOptionalArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("tumor-type", "Optional: Set tumor tissue type to something other than primary tumor (P), i.e. X . Default: Not set (All)").withRequiredArg();
        parser.accepts("interval-bed", "Optional: Specify the path to interval bed file. Default: parsed from " + this.rsConfigXMLPath).withOptionalArg();
        parser.accepts("ref-fasta", "Optional: Specify the path to reference human genome fasta. Default: " + this.refGenome).withOptionalArg();
        parser.accepts("rsconfig-file", "Optional: specify location of .xml file which should be used to cinfigure references. Default: " + this.rsConfigXMLPath).withOptionalArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList(BAM_METATYPE));
        this.setHeadersToGroupBy(Arrays.asList(Header.FILE_SWA));

        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);

        //Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (this.options.has("template-type")) {
            if (!options.hasArgument("template-type")) {
                Log.error("--template-type requires an argument, EX");
            } else {
                this.templateType = options.valueOf("template-type").toString();
                if (!this.templateType.equals("EX")) {
                    Log.stderr("NOTE THAT ONLY EX template-type SUPPORTED, WE CANNOT GUARANTEE MEANINGFUL RESULTS WITH OTHER TEMPLATE TYPES");
                    rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                }
            }
        }
        if (this.options.has("ref-fasta")) {
            this.refGenome = options.valueOf("ref-fasta").toString();
        }

        if (this.options.has("interval-bed")) {
            this.intervalBed = options.valueOf("interval-bed").toString();
        }
        if (this.options.has("rsconfig-file")) {
            this.rsConfigXMLPath = options.valueOf("rsconfig-file").toString();
        }

        try {
            rs = new Rsconfig(new File(this.rsConfigXMLPath));
        } catch (Exception e) {
            Log.error("Rsconfg file did not load properly, exception stack trace:\n" + e.getStackTrace());
            rv.setExitStatus(ReturnValue.FAILURE);
            return rv;
        }

        return rv;
    }

    /**
     * Final check
     *
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     *
     * @return
     */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        /**
         * GP-1627: Update to handle multiple tumour bams matching to one normal
         *
         */
        boolean haveNorm = false;
        boolean haveTumr = false;
        String extTum = null;
        List<String> match = new ArrayList<String>(); // array list to contain one tumour normal pair
        for (String p : filePaths) {
            if (null != this.duplicates && this.duplicates.contains(p)) {
                Log.stderr("File [" + p + "] has a name that cannot be disambiguated in current set, will skip it");
                continue;
            }
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                String tt = bs.tissueType;
                if (!tt.isEmpty() && !tt.equals("R")) {
                    extTum = bs.getRootSampleName();
                    this.tumourFilePath = p;
                    match.add(p); // GP-1627:first element of match array list contains tumour bam file path
                    haveTumr = true;
                    break;
                }
            }
        }

        // Check for duplicate file names and exclude them from analysis
        this.duplicates = detectDuplicates(commaSeparatedFilePaths);

        for (String p : filePaths) {
            if (null != this.duplicates && this.duplicates.contains(p)) {
                Log.stderr("File [" + p + "] has a name that cannot be disambiguated in current set, will skip it");
                continue;
            }
            if (p != this.tumourFilePath) { // check if this path is not tumour bam

                for (BeSmall bs : fileSwaToSmall.values()) {
                    if (!bs.getPath().equals(p)) {
                        continue;
                    }
                    if (bs.getPath().equals(this.tumourFilePath)) { // check if this path is not tumour bam
                        continue;
                    }

                    String rootSampleName = bs.getRootSampleName();
                    String tt = bs.getTissueType();
                    // GP-1627: update match list with normal bam if, root-sample-name for normal == tumour 
                    if (!tt.isEmpty() && tt.equals("R") && rootSampleName.equals(extTum)) {
                        this.normalFilePath = p;
                        match.add(this.normalFilePath);
                        haveNorm = true;
                    }
                }
            } else {
//                Log.debug("move to next");
                continue;
            }
        }
        // size of match list should be exactly == 2
        if (match.size() != 2) {
            Log.error("Does not contain a Tumour Normal pair " + match.size());
        }

        if (haveNorm && haveTumr) {
            // GP-1627: check if match contains a tumour-normal pair
            List<String> accessions = Arrays.asList(commaSeparatedParentAccessions.split(","));
            List<String> files = Arrays.asList(filePaths);
            Map<String, String> fileAccessionMap = new HashMap<String, String>();
            Iterator<String> i1 = files.iterator();
            Iterator<String> i2 = accessions.iterator();
            while (i1.hasNext() && i2.hasNext()) {
                fileAccessionMap.put(i1.next(), i2.next());
            }
            Iterator<String> m = match.iterator();
            List<String> newCommaSeparatedFileArray = new ArrayList<String>();
            List<String> newCommaSeparatedParentAccessionArray = new ArrayList<String>();
            while (m.hasNext()) {
                int index = files.indexOf(m.next());
                newCommaSeparatedFileArray.add(files.get(index));
                newCommaSeparatedParentAccessionArray.add(accessions.get(index));
            }
            // GP-1627: update commaSeparatedFilePaths and commaSeparatedParentAccessions
            this.commaSeparatedFilePaths = String.join(",", newCommaSeparatedFileArray);
            this.commaSeparatedParentAccessions = String.join(",", newCommaSeparatedParentAccessionArray);
            // GP-1627: return the updated tumour-normal filepaths and accessions
            return super.doFinalCheck(this.commaSeparatedFilePaths, this.commaSeparatedParentAccessions);
        }
        String absent = haveNorm ? "Tumor" : "Normal";
        Log.error("Data for " + absent + " tissue are not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
//        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
        String resequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
        String target_bed = null;

        if (null == currentTissueType) {
            return false; // we need only those which have their tissue type set
        }

        // Filter the data of a different template type if filter is specified
        if (!this.templateType.equalsIgnoreCase(currentTtype)) {
            Log.warn("Excluding file with SWID = [" + returnValue.getAttribute(Header.FILE_SWA.getTitle())
                    + "] due to template type/geo_library_source_template_type = [" + currentTtype + "]");
            return false;
        }

        // Do not process tumor tissues of type that doesn't match set parameter
        if (null != this.tumorType) {
            if (!currentTissueType.equals("R") && !currentTissueType.equals(this.tumorType)) {
                return false;
            }
        }

//        Log.debug(rs.get(resequencingType, currentTtype, "interval_file"));
        try {
            target_bed = rs.get(currentTtype, resequencingType, "interval_file").toString();
            if (target_bed != null) {
                this.intervalBed = target_bed;
            } else {
                Log.info("No interval file found for this run; Please re-try with --interval-bed <path to interval file>");
            }
        } catch (Exception NullPointerException) {
            Log.info("No interval file found for this run; Please re-try with --interval-bed <path to interval file>");
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Log.debug("Number of files from file provenance = " + vals.size());

        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(BAM_METATYPE)) {
                        metatypeOK = true;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }

            if (!metatypeOK) {
                continue; // Go to the next value
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();
        /**
         * *
         * GP-1627: Updated to handle multiple tumour files matching to single
         * normal bam Logic: 1. iterate over all tumour bams: 2. For each tumour
         * bam create a list with [tumour bam, normal1, normal2, normal3 .. ]
         *
         *
         */
        // get allTumourOnlyFiles
        Map<String, List<ReturnValue>> tnMap = getTumorNormalMap(newValues);
        return tnMap;
    }
    
    protected List<ReturnValue> getTumourRVs(List<ReturnValue> newValues){
        List<ReturnValue> tumourMap = new ArrayList<ReturnValue> ();
        for (ReturnValue rV : newValues){
//            String sampleName = rV.getAttribute(Header.SAMPLE_NAME.getTitle());
//            ArrayList<FileMetadata> fileMetadata = rV.getFiles();
//            String currentTemplatetype = rV.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            String currentTissueType = rV.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
//            String rootSampleName = rV.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
            if (!currentTissueType.equals("R")){
                tumourMap.add(rV);
            }
        }
        return tumourMap;
    }
    
    
    protected List<ReturnValue> getNormalRVs(List<ReturnValue> newValues){
        List<ReturnValue> normalMap = new ArrayList<ReturnValue> ();
        for (ReturnValue rV : newValues){
            String currentTissueType = rV.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_tissue_type");
            if (currentTissueType.equals("R")){
                normalMap.add(rV);
            }
        }
        return normalMap;
    }
    
    protected Map<String, List<ReturnValue>> getTumorNormalMap(List<ReturnValue> newValues){
        Map<String, List<ReturnValue>> tnMap = new HashMap<String, List<ReturnValue>> ();
        List<ReturnValue> normalBamAttr = getNormalRVs(newValues);
        List<ReturnValue> tumorBamAttr = getTumourRVs(newValues);
        for (ReturnValue tRV : tumorBamAttr){
            // define sample name
            String tumourSampleName = tRV.getAttribute(Header.SAMPLE_NAME.getTitle());
            String tumourRootSampleName = tRV.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
            for (ReturnValue nRV : normalBamAttr){
                String normalSampleName = nRV.getAttribute(Header.SAMPLE_NAME.getTitle());
                String normalRootSampleName = nRV.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
                if (normalRootSampleName.equals(tumourRootSampleName)){
                    String mapKey = normalSampleName+"_"+tumourSampleName;
                    List<ReturnValue> tnList = new ArrayList<ReturnValue> ();
                    tnList.add(tRV);
                    tnList.add(nRV);
                    tnMap.put(mapKey, tnList);
                }
            }
        }
        return tnMap;
        
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String inputNormFile = this.normalFilePath;
        String inputTumrFile = this.tumourFilePath;
        this.sampleName = FilenameUtils.getBaseName(inputTumrFile);
//        String intervalBed = this.rs.get(this.templateType, templateType, queue)

        Map<String, String> iniFileMap = super.modifyIniFile(this.commaSeparatedFilePaths, this.commaSeparatedParentAccessions);

        iniFileMap.put("input_files_normal", inputNormFile);
        iniFileMap.put("input_files_tumor", inputTumrFile);
        iniFileMap.put("data_dir", "data");
        iniFileMap.put("interval_bed", this.intervalBed);
        if (!this.queue.isEmpty()) {
            iniFileMap.put("queue", this.queue);
        }
        iniFileMap.put("external_name", this.sampleName);
        iniFileMap.put("ref_fasta", this.refGenome);
        return iniFileMap;
    }

    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(SequenzaDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String tubeID = null;

        private String extName = null;
        private String groupID = null;
        private String groupDescription = null;
        private String rootSampleName = null;

        public String getRootSampleName() {
            return rootSampleName;
        }

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            extName = rv.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_external_name");
            rootSampleName = rv.getAttribute(Header.ROOT_SAMPLE_NAME.getTitle());
            //fa.getLimsValue(Lims.TUBE_ID);
            if (null == extName || extName.isEmpty()) {
                extName = rootSampleName;
            }
            groupID = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == groupDescription || groupDescription.isEmpty()) {
                groupDescription = "NA";
            }
            StringBuilder gba = new StringBuilder(fa.getDonor());
            gba.append(":").append(fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));

            String trs = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
            if (null != trs && !trs.isEmpty()) {
                gba.append(":").append(trs);
            }

            groupByAttribute = gba.toString() + ":" + extName + ":" + groupID; // grouping issue sequenza decider; generates correct ini file but lists too many files
            path = rv.getFiles().get(0).getFilePath() + "";

        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getTissueType() {
            return tissueType;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public String getTubeID() {
            return tubeID;
        }

        public String getExtName() {
            return extName;
        }

        public String getGroupID() {
            return groupID;
        }

        public String getGroupDescription() {
            return groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }

    public static List<String> detectDuplicates(String commaSeparatedFilePaths) {

        String[] filePaths = commaSeparatedFilePaths.split(",");
        List<String> list = new ArrayList<String>();
        List<String> checker = new ArrayList<String>();

        for (String path : filePaths) {
            String baseName = makeBasename(path, ".bam");

            if (checker.contains(baseName) && !list.contains(path)) {
                list.add(path);
            } else {
                checker.add(baseName);
            }
        }

        return list.isEmpty() ? null : list;

    }

    /**
     * Utility function
     *
     * @param path
     * @param extension
     *
     * @return
     */
    public static String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }
}
