package org.renci.canvas.primer.refseq.commands;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.UUID;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.collections.CollectionUtils;
import org.renci.canvas.primer.commons.FTPFactory;
import org.renci.gbff.GBFFFilter;
import org.renci.gbff.GBFFManager;
import org.renci.gbff.filter.GBFFAndFilter;
import org.renci.gbff.filter.GBFFFeatureSourceOrganismNameFilter;
import org.renci.gbff.filter.GBFFFeatureTypeNameFilter;
import org.renci.gbff.filter.GBFFSequenceAccessionPrefixFilter;
import org.renci.gbff.filter.GBFFSourceOrganismNameFilter;
import org.renci.gbff.model.Sequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DownloadRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(DownloadRunnable.class);

    private static final GBFFManager gbffMgr = GBFFManager.getInstance();

    public DownloadRunnable() {
        super();
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {

            Path outputPath = Paths.get(System.getProperty("karaf.data"), "refseq");
            File refseqDir = outputPath.toFile();

            File gbffDir = new File(refseqDir, "gbff");
            gbffDir.mkdirs();

            List<File> gbffFiles = FTPFactory.ncbiDownloadFiles(gbffDir, "/refseq/release/vertebrate_mammalian", "vertebrate_mammalian",
                    "rna.gbff.gz");

            File alignmentDir = new File(refseqDir, "alignment");
            alignmentDir.mkdirs();

            FTPFactory.ncbiDownloadFiles(alignmentDir, "/refseq/H_sapiens/alignments", "GCF", "refseq_alignments.gff3");

            List<GBFFFilter> filters = Arrays
                    .asList(new GBFFFilter[] { new GBFFSequenceAccessionPrefixFilter(Arrays.asList(new String[] { "NM_", "XM_" })),
                            new GBFFSourceOrganismNameFilter("Homo sapiens"), new GBFFFeatureSourceOrganismNameFilter("Homo sapiens"),
                            new GBFFFeatureTypeNameFilter("CDS"), new GBFFFeatureTypeNameFilter("source") });

            GBFFAndFilter gbffFilter = new GBFFAndFilter(filters);

            List<Sequence> sequences = new ArrayList<>();

            for (File f : gbffFiles) {

                logger.info("parsing GenBankFlatFile: {}", f.getAbsolutePath());
                List<Sequence> sequenceList = gbffMgr.deserialize(gbffFilter, false, f);

                if (CollectionUtils.isNotEmpty(sequenceList)) {
                    logger.info("sequences found: {}", sequenceList.size());
                    for (Sequence sequence : sequenceList) {
                        sequences.add(sequence);
                        if ((sequences.size() % 500) == 0) {

                            File tmpFile = new File(refseqDir, String.format("%s.ser.gz", UUID.randomUUID().toString()));
                            try (FileOutputStream fos = new FileOutputStream(tmpFile);
                                    GZIPOutputStream gzipos = new GZIPOutputStream(fos, Double.valueOf(Math.pow(2, 14)).intValue());
                                    ObjectOutputStream oos = new ObjectOutputStream(gzipos)) {
                                logger.info("serializing: {}", tmpFile.getAbsolutePath());
                                oos.writeObject(sequences);
                            }
                            sequences.clear();

                        }

                    }

                }

                f.delete();
            }

            File tmpFile = new File(refseqDir, String.format("%s.ser.gz", UUID.randomUUID().toString()));
            try (FileOutputStream fos = new FileOutputStream(tmpFile);
                    GZIPOutputStream gzipos = new GZIPOutputStream(fos, Double.valueOf(Math.pow(2, 14)).intValue());
                    ObjectOutputStream oos = new ObjectOutputStream(gzipos)) {
                oos.writeObject(sequences);
            }
            sequences.clear();

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }

    }

}
