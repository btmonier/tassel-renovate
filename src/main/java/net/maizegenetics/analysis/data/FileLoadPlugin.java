/*
 * FileLoadPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.analysis.gobii.GOBIIPlugin;
import net.maizegenetics.dna.map.TOPMUtils;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.ReadSequenceAlignmentUtils;
import net.maizegenetics.dna.snp.io.*;
import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.*;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.taxa.distance.DistanceMatrixUtils;
import net.maizegenetics.taxa.distance.ReadDistanceMatrix;
import net.maizegenetics.taxa.tree.NewickUtils;
import net.maizegenetics.util.HDF5TableReport;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.TableReportUtils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

/**
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class FileLoadPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(FileLoadPlugin.class);

    private PluginParameter<TasselFileType> myFileType = new PluginParameter.Builder<>("format", TasselFileType.Unknown, TasselFileType.class)
            .description("Import file format")
            .objectListSingleSelect()
            .range(TasselFileType.values())
            .build();

    private PluginParameter<Boolean> mySortPositions = new PluginParameter.Builder<>("sortPositions", false, Boolean.class)
            .description("Whether to sort genotype positions if that's possible.")
            .dependentOnParameter(myFileType, new Object[]{TasselFileType.Unknown, TasselFileType.Hapmap, TasselFileType.HapmapDiploid, TasselFileType.VCF, TasselFileType.Plink})
            .build();

    private PluginParameter<Boolean> myKeepDepth = new PluginParameter.Builder<>("keepDepth", false, Boolean.class)
            .description("Whether to keep depth if that's possible.")
            .dependentOnParameter(myFileType, new Object[]{TasselFileType.Unknown, TasselFileType.VCF})
            .build();

    private String[] myOpenFiles = null;
    private PlinkLoadPlugin myPlinkLoadPlugin = null;
    private ProjectionLoadPlugin myProjectionLoadPlugin = null;
    private ProjectPcsAndRunModelSelectionPlugin myProjectPcsAndRunModelSelectionPlugin = null;
    private GOBIIPlugin myGOBIIPlugin = null;
    private final JFileChooser myOpenFileChooser;
    private final boolean myHeadless;

    public enum TasselFileType {

        SqrMatrix("Square Matrix"), Sequence("Sequence"), Unknown("Make Best Guess"),
        Fasta("Fasta"), Hapmap("Hapmap"), HapmapLIX("Hapmap LIX"),
        Plink("Plink"), Phenotype("Phenotype"), ProjectionAlignment("Projection Genotype"),
        ProjectPCsandRunModelSelection("Project PCs"),
        Phylip_Seq("Phylip (Sequential)"), Phylip_Inter("Phylip (Interleaved)"), Table("Table"),
        Serial("Serial"), HapmapDiploid("Hapmap Diploid"), Newick("Newick"), VCF("VCF"),
        HDF5("HDF5"), TOPM("TOPM"), HDF5Schema("HDF5 Schema"), Filter("Filter"),
        NumericGenotype("Numeric Genotype"), TaxaList("Taxa List"), PositionList("Position List"),
        SqrMatrixRaw("Raw MultiBLUP Matrix"), SqrMatrixBin("Binary MultiBLUP Matrix"),
        GOBII("GOBII"), Depth("Depth"), ReferenceProbability("Numeric Genotype"), Report("Report"),
        PlinkPhenotype("Plink Phenotype"), SqrMatrixDARwinDIS("DARwin DIS"), Avro("Avro"),
        Flapjack("Flapjack");

        private final String myText;

        TasselFileType(String text) {
            myText = text;
        }

        @Override
        public String toString() {
            return myText;
        }
    }

    public static final String FILE_EXT_HAPMAP = ".hmp.txt";
    public static final String FILE_EXT_HAPMAP_GZ = ".hmp.txt.gz";
    public static final String FILE_EXT_HAPMAP_GZ_LIX = FILE_EXT_HAPMAP_GZ + LineIndexBuilder.LINE_INDEX_FILE_EXTENSION;
    public static final String FILE_EXT_PLINK_MAP = ".plk.map";
    public static final String FILE_EXT_PLINK_PED = ".plk.ped";
    public static final String FILE_EXT_SERIAL_GZ = ".serial.gz";
    public static final String FILE_EXT_HDF5 = ".h5";
    public static final String FILE_EXT_VCF = ".vcf";
    public static final String FILE_EXT_TOPM = ".topm";
    public static final String FILE_EXT_TOPM_H5 = ".topm.h5";
    public static final String FILE_EXT_TOPM_BIN = ".topm.bin";
    public static final String FILE_EXT_TOPM_TEXT = ".topm.txt";
    public static final String FILE_EXT_FASTA = ".fasta";
    public static final String FILE_EXT_PHYLIP = ".phy";
    public static final String FILE_EXT_NEWICK = ".nwk";

    /**
     * Creates a new instance of FileLoadPlugin. This only used by TASSEL GUI to
     * bypass dialog and go straight to file browser. Bypassing the dialog
     * causes it to bypass adding to Data Tree. This constructor tells
     * FileLoadPlugin to add it to the Data Tree Manually.
     */
    public FileLoadPlugin(Frame parentFrame, boolean isInteractive, boolean headless) {
        super(parentFrame, isInteractive);
        if (isInteractive) {
            myOpenFileChooser = new JFileChooser(TasselPrefs.getOpenDir());
            myOpenFileChooser.setMultiSelectionEnabled(true);
        } else {
            myOpenFileChooser = null;
        }
        myHeadless = headless;
    }

    /**
     * Creates a new instance of FileLoadPlugin.
     */
    public FileLoadPlugin(Frame parentFrame, boolean isInteractive) {
        this(parentFrame, isInteractive, false);
    }

    public static Object runPlugin(String filename) {
        return runPluginDataSet(filename).getData(0).getData();
    }

    public static DataSet runPluginDataSet(String filename) {
        FileLoadPlugin flp = new FileLoadPlugin(null, false);
        flp.setTheFileType(TasselFileType.Unknown);
        flp.setOpenFiles(filename);
        return flp.performFunction(null);
    }

    @Override
    protected void preProcessParameters(DataSet input) {

        List<TasselFileType> temp = new ArrayList<>();
        temp.addAll(Arrays.asList(new FileLoadPlugin.TasselFileType[]{
                TasselFileType.Unknown,
                TasselFileType.Hapmap,
                TasselFileType.VCF,
                TasselFileType.Flapjack,
                TasselFileType.Plink,
                TasselFileType.ProjectionAlignment,
                TasselFileType.Sequence,
                TasselFileType.Fasta,
                TasselFileType.Phenotype,
                TasselFileType.SqrMatrix,
                TasselFileType.Table,
                TasselFileType.Newick,
                TasselFileType.TOPM,
                TasselFileType.HDF5,
                TasselFileType.HDF5Schema}));
        myFileType = new PluginParameter<>(myFileType, temp);

        if (!isInteractive() && myFileType.isEmpty() && myFileType.hasPossibleValues()) {
            fileType(TasselFileType.Unknown);
        }

    }

    @Override
    public DataSet processData(DataSet input) {

        myWasCancelled = true;

        if (isInteractive()) {

            if (fileType() == TasselFileType.Plink) {
                if (myPlinkLoadPlugin == null) {
                    myPlinkLoadPlugin = new PlinkLoadPlugin(getParentFrame(), isInteractive());
                    for (PluginListener current : getListeners()) {
                        myPlinkLoadPlugin.addListener(current);
                    }
                }
                myPlinkLoadPlugin.sortPositions(sortPositions());
                return myPlinkLoadPlugin.performFunction(null);
            }

            if (fileType() == TasselFileType.ProjectionAlignment) {
                if (myProjectionLoadPlugin == null) {
                    myProjectionLoadPlugin = new ProjectionLoadPlugin(getParentFrame(), isInteractive());
                    for (PluginListener current : getListeners()) {
                        myProjectionLoadPlugin.addListener(current);
                    }
                }
                return myProjectionLoadPlugin.performFunction(input);
            }

            if (fileType() == TasselFileType.ProjectPCsandRunModelSelection) {
                if (myProjectPcsAndRunModelSelectionPlugin == null) {
                    myProjectPcsAndRunModelSelectionPlugin = new ProjectPcsAndRunModelSelectionPlugin(getParentFrame(), isInteractive());
                    for (PluginListener current : getListeners()) {
                        myProjectPcsAndRunModelSelectionPlugin.addListener(current);
                    }
                }
                return myProjectPcsAndRunModelSelectionPlugin.performFunction(input);
            }

            if (fileType() == TasselFileType.GOBII) {
                if (myGOBIIPlugin == null) {
                    myGOBIIPlugin = new GOBIIPlugin(getParentFrame(), isInteractive());
                    for (PluginListener current : getListeners()) {
                        myGOBIIPlugin.addListener(current);
                    }
                }
                return myGOBIIPlugin.performFunction(input);
            }

            setOpenFiles(getOpenFilesByChooser());

        }

        if ((myOpenFiles == null) || (myOpenFiles.length == 0)) {
            return null;
        }

        List<DataSet> result = new ArrayList<>();
        ArrayList<String> alreadyLoaded = new ArrayList<>();
        for (int i = 0; i < myOpenFiles.length; i++) {

            if (alreadyLoaded.contains(myOpenFiles[i])) {
                continue;
            }

            LocalDateTime time = LocalDateTime.now();
            String timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Start Loading File: " + myOpenFiles[i] + " time: " + timeStr);

            DataSet tds = null;

            if (fileType() == TasselFileType.Unknown) {
                if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ)) {
                    String theIndex = myOpenFiles[i].replaceFirst(FILE_EXT_HAPMAP_GZ, FILE_EXT_HAPMAP_GZ_LIX);
                    if (new File(theIndex).isFile()) {
                        myLogger.info("guessAtUnknowns: type: " + TasselFileType.HapmapLIX);
                        alreadyLoaded.add(myOpenFiles[i]);
                        alreadyLoaded.add(theIndex);
                        GenotypeTable hapmap = BuilderFromHapMapLIX.build(myOpenFiles[i], theIndex);
                        tds = new DataSet(new Datum(Utils.getFilename(myOpenFiles[i], FileLoadPlugin.FILE_EXT_HAPMAP_GZ), hapmap, null), this);
                    } else {
                        myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                        alreadyLoaded.add(myOpenFiles[i]);
                        tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                    }
                } else if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP_GZ_LIX)) {
                    String theHapmap = myOpenFiles[i].replaceFirst(FILE_EXT_HAPMAP_GZ_LIX, FILE_EXT_HAPMAP_GZ);
                    if (new File(theHapmap).isFile()) {
                        myLogger.info("guessAtUnknowns: type: " + TasselFileType.HapmapLIX);
                        alreadyLoaded.add(myOpenFiles[i]);
                        alreadyLoaded.add(theHapmap);
                        GenotypeTable hapmap = BuilderFromHapMapLIX.build(theHapmap, myOpenFiles[i]);
                        tds = new DataSet(new Datum(Utils.getFilename(theHapmap, FileLoadPlugin.FILE_EXT_HAPMAP_GZ), hapmap, null), this);
                    } else {
                        throw new IllegalStateException("Can't find genotype file for index: " + myOpenFiles[i]);
                    }
                } else if (myOpenFiles[i].endsWith(FILE_EXT_HAPMAP)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Hapmap);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.Hapmap);
                } else if ((myOpenFiles[i].endsWith(FILE_EXT_TOPM_H5)) || (myOpenFiles[i].endsWith(FILE_EXT_TOPM))
                        || (myOpenFiles[i].endsWith(FILE_EXT_TOPM_BIN)) || (myOpenFiles[i].endsWith(FILE_EXT_TOPM_TEXT))) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.TOPM);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.TOPM);
                } else if ((myOpenFiles[i].endsWith(".grm.N.bin")) || (myOpenFiles[i].endsWith(".grm.bin"))
                        || (myOpenFiles[i].endsWith(".grm.id"))) {
                    String[] grmFiles = DistanceMatrixUtils.getGRMFilenames(myOpenFiles[i]);
                    if (new File(grmFiles[0]).isFile() && new File(grmFiles[1]).isFile() && new File(grmFiles[2]).isFile()) {
                        myLogger.info("guessAtUnknowns: type: " + TasselFileType.SqrMatrixBin);
                        alreadyLoaded.add(grmFiles[0]);
                        alreadyLoaded.add(grmFiles[1]);
                        alreadyLoaded.add(grmFiles[2]);
                        tds = processDatum(myOpenFiles[i], TasselFileType.SqrMatrixBin);
                    } else if (myOpenFiles[i].endsWith(".grm.N.bin") && new File(grmFiles[4]).isFile() && new File(myOpenFiles[i]).isFile()) {
                        myLogger.info("guessAtUnknowns: type: " + TasselFileType.SqrMatrix);
                        alreadyLoaded.add(myOpenFiles[i]);
                        alreadyLoaded.add(grmFiles[4]);
                        tds = processDatum(grmFiles[4], TasselFileType.SqrMatrix);
                    }
                } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_PED) || myOpenFiles[i].endsWith(FILE_EXT_PLINK_PED + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                    String theMapFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_PED, FILE_EXT_PLINK_MAP);
                    alreadyLoaded.add(myOpenFiles[i]);
                    alreadyLoaded.add(theMapFile);
                    GenotypeTable plink = ImportUtils.readFromPLink(myOpenFiles[i], theMapFile, this, sortPositions());
                    tds = new DataSet(new Datum(Utils.getFilename(myOpenFiles[i], FileLoadPlugin.FILE_EXT_PLINK_PED), plink, null), this);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_PLINK_MAP) || myOpenFiles[i].endsWith(FILE_EXT_PLINK_MAP + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Plink);
                    String thePedFile = myOpenFiles[i].replaceFirst(FILE_EXT_PLINK_MAP, FILE_EXT_PLINK_PED);
                    alreadyLoaded.add(myOpenFiles[i]);
                    alreadyLoaded.add(thePedFile);
                    GenotypeTable plink = ImportUtils.readFromPLink(thePedFile, myOpenFiles[i], this, sortPositions());
                    tds = new DataSet(new Datum(Utils.getFilename(thePedFile, FileLoadPlugin.FILE_EXT_PLINK_PED), plink, null), this);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_SERIAL_GZ)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Serial);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.Serial);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_HDF5)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.HDF5);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.HDF5);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_VCF) || myOpenFiles[i].endsWith(FILE_EXT_VCF + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.VCF);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.VCF);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_PHYLIP) || myOpenFiles[i].endsWith(FILE_EXT_PHYLIP + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Sequence);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.Sequence);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_FASTA) || myOpenFiles[i].endsWith(FILE_EXT_FASTA + ".gz")) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Fasta);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.Fasta);
                } else if (myOpenFiles[i].endsWith(FILE_EXT_NEWICK)) {
                    myLogger.info("guessAtUnknowns: type: " + TasselFileType.Newick);
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = processDatum(myOpenFiles[i], TasselFileType.Newick);
                } else {
                    alreadyLoaded.add(myOpenFiles[i]);
                    tds = guessAtUnknowns(myOpenFiles[i]);
                }
            } else {
                alreadyLoaded.add(myOpenFiles[i]);
                tds = processDatum(myOpenFiles[i], fileType());
            }

            time = LocalDateTime.now();
            timeStr = time.format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            if (tds != null) {
                myLogger.info("Finished Loading File: " + myOpenFiles[i] + " time: " + timeStr);
                GenotypeSummaryPlugin.printSimpleSummary(tds);
                myWasCancelled = false;
                result.add(tds);
                if (myHeadless) {
                    fireDataSetReturned(new PluginEvent(tds, FileLoadPlugin.class));
                }
            } else {
                myLogger.info("Nothing Loaded for File: " + myOpenFiles[i] + " time: " + timeStr);
            }

        }

        return DataSet.getDataSet(result, this);

    }

    private DataSet guessAtUnknowns(String filename) {

        TasselFileType guess = TasselFileType.Table;

        try (BufferedReader br = Utils.getBufferedReader(filename)) {

            String line1 = br.readLine();
            while (line1 != null) {
                line1 = line1.trim();
                if (!line1.isEmpty()) {
                    break;
                }
                line1 = br.readLine();
            }
            if (line1 == null) {
                throw new IllegalArgumentException("FileLoadPlugin: guessAtUnknowns: File is empty: " + filename);
            }
            String[] sval1 = line1.split("\\s");
            String line2 = br.readLine().trim();
            String[] sval2 = line2.split("\\s");
            if (line1.startsWith("{")) {
                String temp;
                if (sval1.length > 1) {
                    temp = sval1[1];
                } else {
                    temp = line2;
                }
                if (temp.startsWith("\"TaxaList\"")) {
                    guess = TasselFileType.TaxaList;
                } else if (temp.startsWith("\"PositionList\"")) {
                    guess = TasselFileType.PositionList;
                } else if (temp.startsWith("\"Filter\"")) {
                    guess = TasselFileType.Filter;
                }
            } else if (line1.startsWith("##")) {
                String matrixStr = "##" + DistanceMatrixBuilder.MATRIX_TYPE;
                if (line1.startsWith(matrixStr) || line2.startsWith(matrixStr)) {
                    guess = TasselFileType.SqrMatrix;
                } else {
                    String line = br.readLine();
                    while ((line != null) && (line.startsWith("##"))) {
                        if (line.startsWith(matrixStr)) {
                            guess = TasselFileType.SqrMatrix;
                            break;
                        }
                        line = br.readLine();
                    }
                }
            } else if (line1.startsWith("<") || line1.startsWith("#")) {
                boolean isTrait = false;
                boolean isMarker = false;
                boolean isNumeric = false;
                Pattern tagPattern = Pattern.compile("[<>\\s]+");
                String[] info1 = tagPattern.split(line1);
                String[] info2 = tagPattern.split(line2);
                if (info1.length > 1) {
                    if (info1[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info1[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info1[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    } else if (info1[1].toUpperCase().startsWith("PHENO")) {
                        isTrait = true;
                    }
                }
                if (info2.length > 1) {
                    if (info2[1].toUpperCase().startsWith("MARKER")) {
                        isMarker = true;
                    } else if (info2[1].toUpperCase().startsWith("TRAIT")) {
                        isTrait = true;
                    } else if (info2[1].toUpperCase().startsWith("NUMER")) {
                        isNumeric = true;
                    }
                } else {
                    guess = null;
                    String inline = br.readLine();
                    while (guess == null && inline != null && (inline.startsWith("#") || inline.startsWith("<"))) {
                        if (inline.startsWith("<")) {
                            String[] info = tagPattern.split(inline);
                            if (info[1].toUpperCase().startsWith("MARKER")) {
                                isMarker = true;
                            } else if (info[1].toUpperCase().startsWith("TRAIT")) {
                                isTrait = true;
                            } else if (info[1].toUpperCase().startsWith("NUMER")) {
                                isNumeric = true;
                            }
                        }
                    }
                }
                if (isTrait) {
                    guess = TasselFileType.Phenotype;
                } else if (isMarker && isNumeric) {
                    guess = TasselFileType.NumericGenotype;
                } else {
                    myLogger.warn("Line1: " + line1);
                    myLogger.warn("Line2: " + line2);
                    throw new IOException("Improperly formatted header. Data will not be imported for file: " + filename);
                }
            } else if ((line1.startsWith(">")) || (line1.startsWith(";"))) {
                guess = TasselFileType.Fasta;
            } else if (sval1.length == 1) {
                guess = TasselFileType.SqrMatrix;
            } else if ((line1.startsWith("#Nexus")) || (line1.startsWith("#NEXUS")) || (line1.startsWith("CLUSTAL"))
                    || ((sval1.length == 2) && (sval2.length == 2))) {
                guess = TasselFileType.Sequence;
            }

            myLogger.info("guessAtUnknowns: type: " + guess);
            return processDatum(filename, guess);

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("FileLoadPlugin: Problem loading file: " + filename + ".  Error: " + e.getMessage());
        }

    }

    private DataSet processDatum(String inFile, TasselFileType theFT) {

        Object result = null;
        String suffix = null;
        try {
            switch (theFT) {
                case Hapmap: {
                    suffix = FILE_EXT_HAPMAP;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_HAPMAP_GZ;
                    }
                    result = ImportUtils.readFromHapmap(inFile, this, sortPositions());
                    break;
                }
                case HDF5: {
                    IHDF5Reader reader = HDF5Factory.openForReading(inFile);
                    boolean t4HDF5 = HDF5Utils.isTASSEL4HDF5Format(HDF5Factory.openForReading(inFile));
                    reader.close();
                    if (t4HDF5) {
                        String newInfile = inFile.replace(".h5", ".t5.h5");
                        if (new File(newInfile).exists()) {
                            String message = "This file is TASSEL 4 HDF5 format. It looks like it has already been converted to TASSEL 5. Using file: " + newInfile;
                            if (isInteractive()) {
                                DialogUtils.showWarning(message, getParentFrame());
                            } else {
                                myLogger.warn(message);
                            }
                        } else {
                            String message = "This file is TASSEL 4 HDF5 format. It will be converted to TASSEL 5 "
                                    + "HDF5 format with name: " + newInfile + ".  This may take a few minutes.";
                            if (isInteractive()) {
                                DialogUtils.showWarning(message, getParentFrame());
                            } else {
                                myLogger.warn(message);
                            }
                            MigrateHDF5FromT4T5.copyGenotypes(inFile, newInfile);
                        }

                        inFile = newInfile;
                    }
                    suffix = FILE_EXT_HDF5;
                    result = ImportUtils.readGuessFormat(inFile);
                    break;
                }
                case HDF5Schema: {
                    suffix = "";
                    result = new HDF5TableReport(inFile);
                    break;
                }
                case VCF: {
                    suffix = FILE_EXT_VCF;
                    if (inFile.endsWith(".gz")) {
                        suffix = FILE_EXT_VCF + ".gz";
                    }
                    result = ImportUtils.readFromVCF(inFile, this, keepDepth(), sortPositions());
                    break;
                }
                case Sequence: {
                    result = ReadSequenceAlignmentUtils.readBasicAlignments(inFile, 40);
                    break;
                }
                case Fasta: {
                    result = ImportUtils.readFasta(inFile);
                    break;
                }
                case SqrMatrix: {
                    result = ReadDistanceMatrix.readDistanceMatrix(inFile);
                    break;
                }
                case SqrMatrixBin: {
                    result = ReadDistanceMatrix.readBinMultiBlupMatrix(inFile);
                    break;
                }
                case Phenotype: {
                    List<Phenotype> phenotypes = new PhenotypeBuilder().fromFile(inFile).build();
                    if (phenotypes.size() != 1) {
                        throw new IllegalStateException("FileLoadPlugin: processDatum: problem loading phenotype file: " + inFile);
                    }
                    result = phenotypes.get(0);
                    break;
                }
                case NumericGenotype: {
                    result = ReadNumericMarkerUtils.readNumericMarkerFile(inFile);
                    break;
                }
                case TaxaList: {
                    result = JSONUtils.importTaxaListFromJSON(inFile);
                    break;
                }
                case PositionList: {
                    result = JSONUtils.importPositionListFromJSON(inFile);
                    break;
                }
                case Table: {
                    result = TableReportUtils.readDelimitedTableReport(inFile, "\t");
                    break;
                }
                case TOPM: {
                    result = TOPMUtils.readTOPM(inFile);
                    break;
                }
                case Filter: {
                    result = FilterJSONUtils.importJSONToFilter(inFile);
                    break;
                }
                case Newick: {
                    result = NewickUtils.read(inFile);
                    break;
                }
                default: {
                    throw new IllegalStateException("Unknown Format: " + theFT + ".\n  Please check file format or select specific format.");
                }
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("Problem loading file: " + inFile + ".\n  Error: " + e.getMessage());
        }

        if (result != null) {
            Datum td = new Datum(Utils.getFilename(inFile, suffix), result, null);
            return new DataSet(td, this);
        }
        return null;

    }

    /**
     * Provides a open filer that remember the last location something was
     * opened from
     */
    private File[] getOpenFilesByChooser() {
        File[] lopenFiles = null;
        myOpenFileChooser.setVisible(true);
        int returnVal = myOpenFileChooser.showOpenDialog(getParentFrame());
        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            lopenFiles = myOpenFileChooser.getSelectedFiles();
            TasselPrefs.putOpenDir(myOpenFileChooser.getCurrentDirectory().getPath());
        }
        return lopenFiles;
    }

    public String[] getOpenFiles() {
        return myOpenFiles;
    }

    public void setOpenFiles(File[] openFiles) {

        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
            return;
        }

        myOpenFiles = new String[openFiles.length];
        for (int i = 0; i < openFiles.length; i++) {
            myOpenFiles[i] = openFiles[i].getPath();
        }

    }

    public void setOpenFiles(String openFile) {
        if ((openFile == null) || openFile.isEmpty()) {
            myOpenFiles = null;
        } else {
            myOpenFiles = new String[]{openFile};
        }
    }

    public void setOpenFiles(String[] openFiles) {
        if ((openFiles == null) || (openFiles.length == 0)) {
            myOpenFiles = null;
        } else {
            myOpenFiles = openFiles;
        }
    }

    /**
     * Export file format (Default format depends on data being exported)
     *
     * @return Format
     */
    public TasselFileType fileType() {
        return myFileType.value();
    }

    /**
     * Set Format. Export file format (Default format depends on data being
     * exported)
     *
     * @param value Format
     *
     * @return this plugin
     */
    public FileLoadPlugin fileType(TasselFileType value) {
        myFileType = new PluginParameter<>(myFileType, value);
        return this;
    }

    public TasselFileType getTheFileType() {
        return fileType();
    }

    public void setTheFileType(TasselFileType theFileType) {
        fileType(theFileType);
    }

    /**
     * Whether to sort genotype positions if that's possible.
     *
     * @return Sort Positions
     */
    public Boolean sortPositions() {
        return mySortPositions.value();
    }

    /**
     * Set Sort Positions. Whether to sort genotype positions if that's
     * possible.
     *
     * @param value Sort Positions
     *
     * @return this plugin
     */
    public FileLoadPlugin sortPositions(Boolean value) {
        mySortPositions = new PluginParameter<>(mySortPositions, value);
        return this;
    }

    /**
     * Whether to keep depth if that's possible.
     *
     * @return Keep Depth
     */
    public Boolean keepDepth() {
        return myKeepDepth.value();
    }

    /**
     * Set Keep Depth. Whether to keep depth if that's possible.
     *
     * @param value Keep Depth
     *
     * @return this plugin
     */
    public FileLoadPlugin keepDepth(Boolean value) {
        myKeepDepth = new PluginParameter<>(myKeepDepth, value);
        return this;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    @Override
    public ImageIcon getIcon() {
        URL imageURL = FileLoadPlugin.class.getResource("/net/maizegenetics/analysis/images/LoadFile.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    @Override
    public String getButtonName() {
        return "Open As...";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    @Override
    public String getToolTipText() {
        return "Open data from filesystem.";
    }
}
