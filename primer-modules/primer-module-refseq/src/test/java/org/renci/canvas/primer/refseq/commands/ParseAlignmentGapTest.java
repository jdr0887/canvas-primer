package org.renci.canvas.primer.refseq.commands;

import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.junit.Test;
import org.renci.canvas.dao.refseq.model.TranscriptMapsExons;
import org.renci.canvas.dao.refseq.model.TranscriptMapsExonsPK;
import org.renci.gff3.GFF3Manager;
import org.renci.gff3.model.GFF3Record;
import org.renci.gff3.model.StrandType;

public class ParseAlignmentGapTest {

    @Test
    public void test() {

        GFF3Manager gff3Mgr = GFF3Manager.getInstance();

        File downloadsDir = new File("/home/jdr0887/Downloads");

        List<File> alignmentFiles = Arrays.asList(new File(downloadsDir,
                "GCF_000001405.25_knownrefseq_alignments.gff3")/*
                                                                * , new File(downloadsDir, "GCF_000001405.25_modelrefseq_alignments.gff3"),
                                                                * new File(downloadsDir, "GCF_000001405.28_knownrefseq_alignments.gff3"),
                                                                * new File(downloadsDir, "GCF_000001405.28_modelrefseq_alignments.gff3"),
                                                                * new File(downloadsDir, "GCF_000001405.33_knownrefseq_alignments.gff3"),
                                                                * new File(downloadsDir, "GCF_000001405.33_modelrefseq_alignments.gff3")
                                                                */);

        Map<String, List<GFF3Record>> alignmentMap = new HashMap<>();

        for (File alignmentFile : alignmentFiles) {
            List<GFF3Record> results = gff3Mgr.deserialize(alignmentFile);
            if (CollectionUtils.isNotEmpty(results)) {
                for (GFF3Record record : results) {
                    if (!record.getSequenceId().startsWith("NC_")) {
                        continue;
                    }
                    alignmentMap.putIfAbsent(record.getSequenceId(), new ArrayList<>());
                    alignmentMap.get(record.getSequenceId()).add(record);
                }
            }
        }

        // NM_173600.2; M6352 I2 M5901
        // NM_005578.4 M9944 D1 M6377
        // NM_001278267.1 M174 I3 M3 I3 M5
        // NM_000294.2
        // NM_001145103.1 1029 5808; M3186 D3 M1594
        // NM_001278267.1 1209 1396; M174 I3 M3 I3 M5

        List<TranscriptMapsExons> exons = new ArrayList<>();

        for (String sequenceId : alignmentMap.keySet()) {

            List<GFF3Record> records = alignmentMap.get(sequenceId);
            List<GFF3Record> filteredRecords = records.stream().filter(a -> a.getAttributes().get("Target").startsWith("NM_001278267.1"))
                    .collect(Collectors.toList());
            if (CollectionUtils.isEmpty(filteredRecords)) {
                continue;
            }

            String genomeReferenceAccession = sequenceId;
            String strand = filteredRecords.stream().map(a -> a.getStrand().getSymbol()).distinct().collect(Collectors.joining());
            String identity = filteredRecords.stream().map(a -> a.getAttributes().get("identity")).distinct().collect(Collectors.joining());

            List<Integer> contigCoordinates = new ArrayList<>();

            int recordIdx = 1;
            for (GFF3Record record : filteredRecords) {

                String gap = record.getAttributes().get("Gap");

                String target = record.getAttributes().get("Target");
                String[] parts = target.split(" ");

                Integer exonStart = Integer.valueOf(parts[1]);
                Integer exonStop = Integer.valueOf(parts[2]);

                if (StringUtils.isEmpty(gap)) {

                    TranscriptMapsExonsPK exonPK = new TranscriptMapsExonsPK(null, recordIdx++);
                    TranscriptMapsExons exon = new TranscriptMapsExons(exonPK);

                    exon.setContigStart(record.getStart());
                    contigCoordinates.add(record.getStart());

                    exon.setContigEnd(record.getEnd());
                    contigCoordinates.add(record.getEnd());

                    exon.setTranscriptStart(exonStart);
                    exon.setTranscriptEnd(exonStop);

                    exons.add(exon);
                } else {

                    if (record.getStrand().equals(StrandType.POSITIVE)) {

                        Integer contigStart = record.getStart();

                        String[] gapTokens = gap.split(" ");
                        for (String gapToken : gapTokens) {
                            String firstChar = gapToken.substring(0, 1);
                            String length = gapToken.substring(1, gapToken.length());

                            CIGARType type = Arrays.asList(CIGARType.values()).stream().filter(a -> a.getName().equals(firstChar))
                                    .findFirst().orElse(null);
                            CIGARElement element = new CIGARElement(type, Integer.valueOf(length));

                            if (element.getType().equals(CIGARType.MATCH)) {

                                TranscriptMapsExonsPK exonPK = new TranscriptMapsExonsPK(null, recordIdx++);
                                TranscriptMapsExons exon = new TranscriptMapsExons(exonPK);

                                exon.setContigStart(contigStart);
                                exon.setContigEnd(contigStart + element.getLength() - 1);

                                exon.setTranscriptStart(exonStart);
                                exon.setTranscriptEnd(exonStart + element.getLength() - 1);
                                exons.add(exon);

                                contigStart = exon.getContigEnd();
                                exonStart = exon.getTranscriptEnd();
                            }

                            if (element.getType().equals(CIGARType.INSERT)) {
                                contigStart++;
                                exonStart = exonStart + element.getLength() + 1;
                            }

                            if (element.getType().equals(CIGARType.DELETION)) {
                                contigStart = contigStart + element.getLength() + 1;
                                exonStart++;
                            }
                        }
                    } else {

                        Integer contigStart = record.getStart();

                        String[] gapTokens = gap.split(" ");
                        for (String gapToken : gapTokens) {
                            String firstChar = gapToken.substring(0, 1);
                            String length = gapToken.substring(1, gapToken.length());
                            CIGARType type = Arrays.asList(CIGARType.values()).stream().filter(a -> a.getName().equals(firstChar))
                                    .findFirst().orElse(null);
                            CIGARElement element = new CIGARElement(type, Integer.valueOf(length));

                            if (element.getType().equals(CIGARType.MATCH)) {

                                TranscriptMapsExonsPK exonPK = new TranscriptMapsExonsPK(null, recordIdx++);
                                TranscriptMapsExons exon = new TranscriptMapsExons(exonPK);

                                exon.setContigStart(contigStart);
                                exon.setContigEnd(contigStart + element.getLength() - 1);

                                exon.setTranscriptStart(exonStart);
                                exon.setTranscriptEnd(exonStart + element.getLength() - 1);
                                exons.add(exon);

                                contigStart = exon.getContigEnd();
                                exonStart = exon.getTranscriptEnd();
                            }

                            if (element.getType().equals(CIGARType.INSERT)) {
                                contigStart++;
                                exonStart = exonStart + element.getLength() + 1;
                            }

                            if (element.getType().equals(CIGARType.DELETION)) {
                                contigStart = contigStart + element.getLength() + 1;
                                exonStart++;
                            }
                        }

                    }

                }

            }

            exons.forEach(a -> System.out.println(a.toString()));

            for (TranscriptMapsExons exon : exons) {
                assertTrue(exon.getContigEnd() - exon.getContigStart() == exon.getTranscriptEnd() - exon.getTranscriptStart());
            }

        }

    }

    @Test
    public void scratch() {

        GFF3Manager gff3Mgr = GFF3Manager.getInstance();

        File downloadsDir = new File("/home/jdr0887/Downloads");

        List<File> alignmentFiles = Arrays.asList(new File(downloadsDir, "GCF_000001405.28_knownrefseq_alignments.gff3"),
                // new File(downloadsDir, "GCF_000001405.28_modelrefseq_alignments.gff3"),
                new File(downloadsDir, "GCF_000001405.33_knownrefseq_alignments.gff3")
        // new File(downloadsDir, "GCF_000001405.33_modelrefseq_alignments.gff3")

        );

        Map<String, List<GFF3Record>> alignmentMap = new HashMap<>();

        for (File alignmentFile : alignmentFiles) {
            List<GFF3Record> results = gff3Mgr.deserialize(alignmentFile);
            if (CollectionUtils.isNotEmpty(results)) {
                for (GFF3Record record : results) {
                    if (!record.getSequenceId().startsWith("NC_")) {
                        continue;
                    }
                    alignmentMap.putIfAbsent(record.getSequenceId(), new ArrayList<>());
                    alignmentMap.get(record.getSequenceId()).add(record);
                }
            }
        }

        try (FileWriter fw = new FileWriter(new File("/tmp", "asdf.sql")); BufferedWriter bw = new BufferedWriter(fw)) {

            for (String sequenceId : alignmentMap.keySet()) {

                List<GFF3Record> records = alignmentMap.get(sequenceId);
                for (GFF3Record record : records) {

                    String gap = record.getAttributes().get("Gap");

                    String target = record.getAttributes().get("Target");
                    String[] parts = target.split(" ");

                    String transcriptVersionId = parts[0];
                    Integer exonStart = Integer.valueOf(parts[1]);
                    Integer exonStop = Integer.valueOf(parts[2]);

                    if (StringUtils.isNotEmpty(gap)) {

                        // bw.write(String.format(
                        // "select * from refseq.transcr_maps_exons a join refseq.transcr_maps b on a.refseq_transcr_maps_id =
                        // b.refseq_transcr_maps_id where b.genome_ref_id = 4 and b.refseq_transcr_ver_id = '%s' and transcr_start = %s and
                        // transcr_end = %s and b.seq_ver_accession = '%s';",
                        // transcriptVersionId, exonStart, exonStop, record.getSequenceId()));
                        bw.write(String.format(
                                "update refseq.transcr_maps_exons a set gap = '%s' from refseq.transcr_maps b where a.refseq_transcr_maps_id = b.refseq_transcr_maps_id and b.genome_ref_id = 4 and b.refseq_transcr_ver_id = '%s' and a.transcr_start = %s and a.transcr_end = %s and b.seq_ver_accession = '%s';",
                                gap, transcriptVersionId, exonStart, exonStop, record.getSequenceId()));
                        bw.newLine();
                        bw.flush();
                    }

                }

            }

        } catch (IOException e) {
            e.printStackTrace();
        }

    }

}
