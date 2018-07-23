package org.renci.canvas.primer.refseq.commands;

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.annotation.model.AnnotationGene;
import org.renci.canvas.dao.annotation.model.AnnotationGeneExternalId;
import org.renci.canvas.dao.annotation.model.AnnotationGeneExternalIdPK;
import org.renci.canvas.dao.annotation.model.AnnotationGeneSynonym;
import org.renci.canvas.dao.annotation.model.AnnotationGeneSynonymPK;
import org.renci.canvas.dao.ref.model.GenomeRef;
import org.renci.canvas.dao.ref.model.GenomeRefSeq;
import org.renci.canvas.dao.refseq.model.CDSECNumber;
import org.renci.canvas.dao.refseq.model.CDSECNumberPK;
import org.renci.canvas.dao.refseq.model.CDSTranslationException;
import org.renci.canvas.dao.refseq.model.CDSTranslationExceptionPK;
import org.renci.canvas.dao.refseq.model.Feature;
import org.renci.canvas.dao.refseq.model.FeatureType;
import org.renci.canvas.dao.refseq.model.GroupingType;
import org.renci.canvas.dao.refseq.model.RefSeqCodingSequence;
import org.renci.canvas.dao.refseq.model.RefSeqGene;
import org.renci.canvas.dao.refseq.model.RegionGroup;
import org.renci.canvas.dao.refseq.model.RegionGroupRegion;
import org.renci.canvas.dao.refseq.model.RegionGroupRegionPK;
import org.renci.canvas.dao.refseq.model.Transcript;
import org.renci.canvas.dao.refseq.model.TranscriptMaps;
import org.renci.canvas.dao.refseq.model.TranscriptMapsExons;
import org.renci.canvas.dao.refseq.model.TranscriptMapsExonsPK;
import org.renci.canvas.dao.refseq.model.TranscriptRefSeqVersion;
import org.renci.canvas.dao.refseq.model.TranscriptRefSeqVersionPK;
import org.renci.canvas.primer.commons.FTPFactory;
import org.renci.canvas.primer.commons.UpdateDiagnosticResultVersionCallable;
import org.renci.gbff.model.Sequence;
import org.renci.gbff.model.TranslationException;
import org.renci.gff3.GFF3Manager;
import org.renci.gff3.model.GFF3Record;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "refseq", name = "persist", description = "Download & persist RefSeq data")
@Service
public class PersistAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(PersistAction.class);

    private static final Pattern locationPattern = Pattern.compile("(?<start>\\d+)\\.+?(?<stop>\\d+)?");

    @Reference
    private CANVASDAOBeanService canvasDAOBeanService;

    public PersistAction() {
        super();
    }

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");

        Executors.newSingleThreadExecutor().execute(() -> {
            long start = System.currentTimeMillis();

            try {

                List<String> remoteFileNames = FTPFactory.ncbiListRemoteFiles("/refseq/release/release-catalog", "release");
                if (CollectionUtils.isEmpty(remoteFileNames)) {
                    logger.error("No remote files found to get version");
                    return;
                }

                Path alignmentPath = Paths.get(System.getProperty("karaf.data"), "refseq", "alignment");
                File alignmentDir = alignmentPath.toFile();
                List<File> alignmentFiles = Arrays.asList(alignmentDir.listFiles());

                GFF3Manager gff3Mgr = GFF3Manager.getInstance();
                Map<String, List<GFF3Record>> alignmentMap = new HashMap<>();
                for (File alignmentFile : alignmentFiles) {
                    logger.info("parsing alignment file: {}", alignmentFile.getName());
                    List<GFF3Record> results = gff3Mgr.deserialize(alignmentFile);
                    if (CollectionUtils.isNotEmpty(results)) {
                        for (GFF3Record record : results) {
                            alignmentMap.putIfAbsent(record.getSequenceId(), new ArrayList<>());
                            alignmentMap.get(record.getSequenceId()).add(record);
                        }
                    }
                }

                String first = remoteFileNames.get(0);
                String refseqVersion = first.substring(7, first.indexOf("."));
                logger.info("refseqVersion = {}", refseqVersion);

                List<GroupingType> allGroupingTypes = canvasDAOBeanService.getGroupingTypeDAO().findAll();

                GroupingType singleGroupingType = allGroupingTypes.stream().filter(a -> a.getId().equals("single")).findAny().get();

                Path outputPath = Paths.get(System.getProperty("karaf.data"), "refseq");
                File refseqDir = outputPath.toFile();

                List<File> serializedSequenceFiles = Arrays.asList(refseqDir.listFiles());

                logger.info("serializedSequenceFiles.size(): {}", serializedSequenceFiles.size());

                for (File serializedSequenceFile : serializedSequenceFiles) {

                    List<Sequence> sequenceList = null;

                    try (FileInputStream fis = new FileInputStream(serializedSequenceFile);
                            GZIPInputStream gzipis = new GZIPInputStream(fis, Double.valueOf(Math.pow(2, 16)).intValue());
                            ObjectInputStream ois = new ObjectInputStream(gzipis)) {
                        logger.info("deserializing: {}", serializedSequenceFile.getAbsolutePath());
                        sequenceList = (List<Sequence>) ois.readObject();
                    } catch (Exception e) {
                        logger.error(e.getMessage(), e);
                    }

                    logger.info("sequenceList.size(): {}", sequenceList.size());

                    List<String> featureExclusionList = Arrays.asList("CDS", "gene", "STS", "variation", "exon", "source", "precursor_RNA");

                    for (Sequence sequence : sequenceList) {
                        try {
                            List<org.renci.gbff.model.Feature> features = sequence.getFeatures();
                            List<org.renci.gbff.model.Feature> filteredFeatures = features.stream()
                                    .filter(a -> !featureExclusionList.contains(a.getType())).collect(Collectors.toList());
                            for (org.renci.gbff.model.Feature gbffFeature : filteredFeatures) {
                                FeatureType featureType = canvasDAOBeanService.getFeatureTypeDAO().findById(gbffFeature.getType());
                                if (featureType == null) {
                                    canvasDAOBeanService.getFeatureTypeDAO().save(new FeatureType(gbffFeature.getType()));
                                }
                            }
                        } catch (CANVASDAOException e) {
                            logger.error(e.getMessage(), e);
                        }
                    }

                    ExecutorService es = Executors.newFixedThreadPool(3);
                    for (Sequence sequence : sequenceList) {
                        es.submit(() -> {

                            logger.info(sequence.toString());
                            try {

                                Transcript transcript = persistTranscript(refseqVersion, sequence);
                                logger.info(transcript.toString());

                                persistGenes(refseqVersion, transcript, singleGroupingType, sequence.getFeatures());

                                persistCodingSequence(refseqVersion, transcript, singleGroupingType, sequence.getFeatures());

                                persistFeatures(refseqVersion, transcript, allGroupingTypes, sequence.getFeatures());

                                persistMappings(refseqVersion, transcript, alignmentMap);

                            } catch (Exception e) {
                                logger.error(e.getMessage(), e);
                            }

                        });
                    }
                    es.shutdown();
                    if (!es.awaitTermination(1L, TimeUnit.DAYS)) {
                        es.shutdownNow();
                    }

                    serializedSequenceFile.delete();

                }

                UpdateDiagnosticResultVersionCallable callable = new UpdateDiagnosticResultVersionCallable(canvasDAOBeanService);
                callable.setNote(String.format("Pulling latest RefSeq: %s", refseqVersion));
                Executors.newSingleThreadExecutor().submit(callable).get();

            } catch (Exception e) {
                logger.error(e.getMessage(), e);
            }
            long end = System.currentTimeMillis();
            logger.info("duration = {}", String.format("%d seconds", (end - start) / 1000));

        });

        return null;
    }

    private void persistMappings(String refseqVersion, Transcript transcript, Map<String, List<GFF3Record>> map) throws CANVASDAOException {
        logger.debug("ENTERING persistMappings(String, Transcript, List<File>)");

        Map<String, List<GFF3Record>> filteredMap = map.entrySet().stream()
                .filter(a -> CollectionUtils.isNotEmpty(a.getValue().stream()
                        .filter(b -> b.getAttributes().get("Target").startsWith(transcript.getId())).collect(Collectors.toSet())))
                .collect(Collectors.toMap(p -> p.getKey(), p -> p.getValue()));

        for (String sequenceId : filteredMap.keySet()) {
            List<GFF3Record> records = map.get(sequenceId);

            String genomeReferenceAccession = sequenceId;
            String strand = records.stream().map(a -> a.getStrand().getSymbol()).distinct().collect(Collectors.joining());
            String identity = records.stream().map(a -> a.getAttributes().get("identity")).distinct().collect(Collectors.joining());

            if ("null".equals(genomeReferenceAccession)) {
                logger.error("Could not get valid GenomeRefSeq for: {}", genomeReferenceAccession);
                return;
            }

            logger.info("genomeReferenceAccession: {}", genomeReferenceAccession);

            if ("null".equals(strand)) {
                logger.error("Could not get valid strand for: {}", genomeReferenceAccession);
                return;
            }

            logger.info("strand: {}", strand);

            GenomeRefSeq genomeRefSeq = canvasDAOBeanService.getGenomeRefSeqDAO().findById(genomeReferenceAccession);
            if (genomeRefSeq == null) {
                logger.error("GenomeRefSeq not found: {}", genomeReferenceAccession);
                return;
            }
            logger.info(genomeRefSeq.toString());

            List<GenomeRef> foundGenomeRefs = canvasDAOBeanService.getGenomeRefDAO()
                    .findByGenomeRefSeqVersionAccession(genomeRefSeq.getId());
            if (CollectionUtils.isEmpty(foundGenomeRefs)) {
                logger.error("GenomeRef not found: {}", genomeReferenceAccession);
                return;
            }
            GenomeRef genomeRef = foundGenomeRefs.get(0);
            logger.info(genomeRef.toString());

            TranscriptMaps transcriptMaps = new TranscriptMaps();

            transcriptMaps.setTranscript(transcript);

            transcriptMaps.setGenomeRef(genomeRef);
            transcriptMaps.setGenomeRefSeq(genomeRefSeq);
            transcriptMaps.setStrand(strand);
            transcriptMaps.setMapCount(canvasDAOBeanService.getTranscriptMapsDAO().findNextMapCount());

            if (StringUtils.isNotEmpty(identity) && !"null".equals(identity)) {
                transcriptMaps.setIdentity(Double.valueOf(identity) * 100D);
            } else {
                identity = records.stream().map(a -> a.getAttributes().get("pct_identity_gap")).distinct().collect(Collectors.joining());
                transcriptMaps.setIdentity(Double.valueOf(identity));
            }

            transcriptMaps.setScore(transcriptMaps.getIdentity());
            transcriptMaps.setExonCount(records.size());

            List<TranscriptMaps> foundTranscriptMaps = canvasDAOBeanService.getTranscriptMapsDAO().findByExample(transcriptMaps);
            if (CollectionUtils.isNotEmpty(foundTranscriptMaps)) {
                transcriptMaps = foundTranscriptMaps.get(0);
            } else {
                transcriptMaps.setId(canvasDAOBeanService.getTranscriptMapsDAO().save(transcriptMaps));
            }

            List<Integer> contigCoordinates = new ArrayList<>();

            int recordIdx = 1;
            for (GFF3Record record : records) {

                String gap = record.getAttributes().get("Gap");

                String target = record.getAttributes().get("Target");
                String[] parts = target.split(" ");

                String exonStart = parts[1];
                String exonStop = parts[2];

                try {
                    TranscriptMapsExonsPK exonPK = new TranscriptMapsExonsPK(transcriptMaps.getId(), recordIdx++);
                    TranscriptMapsExons exon = new TranscriptMapsExons(exonPK);
                    exon.setTranscriptMaps(transcriptMaps);

                    exon.setContigStart(record.getStart());
                    contigCoordinates.add(record.getStart());

                    exon.setContigEnd(record.getEnd());
                    contigCoordinates.add(record.getEnd());

                    if (StringUtils.isNotEmpty(gap)) {
                        exon.setGap(gap);
                    }

                    exon.setTranscriptStart(Integer.valueOf(exonStart));
                    exon.setTranscriptEnd(Integer.valueOf(exonStop));

                    canvasDAOBeanService.getTranscriptMapsExonsDAO().save(exon);
                    transcriptMaps.getExons().add(exon);
                } catch (Exception e) {
                    logger.error(e.getMessage(), e);
                }

            }

            transcriptMaps.setMinContig(Collections.min(contigCoordinates));
            transcriptMaps.setMaxContig(Collections.max(contigCoordinates));

            canvasDAOBeanService.getTranscriptMapsDAO().save(transcriptMaps);
            logger.info(transcriptMaps.toString());

            transcriptMaps.getExons().stream().forEach(a -> logger.info(a.toString()));

        }

    }

    private void persistFeatures(String refseqVersion, Transcript transcript, List<GroupingType> allGroupingTypes,
            List<org.renci.gbff.model.Feature> features) throws CANVASDAOException {
        logger.debug("ENTERING persistFeatures(String, Transcript, List<GroupingType>, List<Feature>)");

        List<String> featureExclusionList = Arrays.asList("CDS", "gene", "STS", "variation", "exon", "source", "precursor_RNA");
        List<org.renci.gbff.model.Feature> filteredFeatures = features.stream().filter(a -> !featureExclusionList.contains(a.getType()))
                .collect(Collectors.toList());

        for (org.renci.gbff.model.Feature gbffFeature : filteredFeatures) {
            logger.info(gbffFeature.toString());

            FeatureType featureType = canvasDAOBeanService.getFeatureTypeDAO().findById(gbffFeature.getType());

            Map<String, String> qualifiers = gbffFeature.getQualifiers();
            String note = qualifiers.get("note");
            String location = gbffFeature.getLocation();

            List<Pair<String, String>> rangeList = new ArrayList<>();

            GroupingType groupingType;
            if (gbffFeature.getLocation().startsWith("join")) {
                groupingType = allGroupingTypes.stream().filter(a -> a.getId().equals("join")).findAny().get();
                Arrays.asList(location.substring(6, location.length() - 1).split(",")).forEach(a -> {
                    Matcher m = locationPattern.matcher(a);
                    m.find();
                    if (m.matches()) {
                        rangeList.add(Pair.of(m.group("start"), m.group("stop")));
                    }
                });
            } else if (gbffFeature.getLocation().startsWith("order")) {
                Arrays.asList(location.substring(7, location.length() - 1).split(",")).forEach(a -> {
                    Matcher m = locationPattern.matcher(a);
                    m.find();
                    if (m.matches()) {
                        rangeList.add(Pair.of(m.group("start"), m.group("stop")));
                    }
                });
                groupingType = allGroupingTypes.stream().filter(a -> a.getId().equals("order")).findAny().get();
            } else {
                Matcher m = locationPattern.matcher(location);
                m.find();
                if (m.matches()) {
                    rangeList.add(Pair.of(m.group("start"), m.group("stop")));
                }
                groupingType = allGroupingTypes.stream().filter(a -> a.getId().equals("single")).findAny().get();
            }
            logger.debug(groupingType.toString());

            RegionGroup regionGroup = new RegionGroup();
            regionGroup.setTranscript(transcript);
            regionGroup.setGroupingType(groupingType);
            regionGroup.setId(canvasDAOBeanService.getRegionGroupDAO().save(regionGroup));
            logger.info(regionGroup.toString());

            for (Pair<String, String> pair : rangeList) {

                String startLocation = pair.getLeft();
                String stopLocation = pair.getRight();

                String startType = "exact";
                String stopType = "exact";

                if (pair.getLeft().startsWith("<")) {
                    startType = "less-than";
                    startLocation = startLocation.substring(1, startLocation.length());
                } else if (pair.getLeft().startsWith(">")) {
                    startType = "greater-than";
                    startLocation = startLocation.substring(1, startLocation.length());
                }

                if (pair.getRight().startsWith("<")) {
                    stopType = "less-than";
                    stopLocation = stopLocation.substring(1, stopLocation.length());
                } else if (pair.getRight().startsWith(">")) {
                    stopType = "greater-than";
                    stopLocation = stopLocation.substring(1, stopLocation.length());
                }

                RegionGroupRegionPK rgrPK = new RegionGroupRegionPK(Integer.valueOf(startLocation), Integer.valueOf(stopLocation),
                        startType, stopType, regionGroup.getId());

                RegionGroupRegion rgr = new RegionGroupRegion(rgrPK);
                rgr.setRegionGroup(regionGroup);
                logger.info(rgr.toString());
                canvasDAOBeanService.getRegionGroupRegionDAO().save(rgr);

            }

            Feature feature = new Feature(refseqVersion, note);
            feature.setType(featureType);
            feature.setRegionGroup(regionGroup);
            feature.setId(canvasDAOBeanService.getFeatureDAO().save(feature));
            logger.info(feature.toString());

        }

    }

    private void persistCodingSequence(String refseqVersion, Transcript transcript, GroupingType singleGroupingType,
            List<org.renci.gbff.model.Feature> features) throws CANVASDAOException {
        logger.debug("ENTERING persistCodingSequence(String, Transcript, List<Feature>)");
        Optional<org.renci.gbff.model.Feature> optionalCodingSequenceFeature = features.stream().filter(a -> "CDS".equals(a.getType()))
                .findAny();
        if (optionalCodingSequenceFeature.isPresent()) {
            org.renci.gbff.model.Feature codingSequenceFeature = optionalCodingSequenceFeature.get();
            logger.info(codingSequenceFeature.toString());

            Map<String, String> qualifiers = codingSequenceFeature.getQualifiers();
            Map<String, String> dbxrefMap = codingSequenceFeature.getDbXRefs();
            List<TranslationException> translationExceptionList = codingSequenceFeature.getTranslationExceptions();

            String note = qualifiers.get("note");
            String codonStart = qualifiers.get("codon_start");
            String product = qualifiers.get("product");
            String proteinId = qualifiers.get("protein_id");
            String translation = qualifiers.get("translation");
            String description = qualifiers.get("desc");
            String ecNumber = qualifiers.get("EC_number");

            RefSeqCodingSequence refseqCodingSequence = new RefSeqCodingSequence(refseqVersion, proteinId, product,
                    Integer.valueOf(codonStart), description, translation, note);

            List<RefSeqCodingSequence> foundRefSeqCodingSequences = canvasDAOBeanService.getRefSeqCodingSequenceDAO()
                    .findByExample(refseqCodingSequence);

            if (CollectionUtils.isEmpty(foundRefSeqCodingSequences)) {
                refseqCodingSequence.setId(canvasDAOBeanService.getRefSeqCodingSequenceDAO().save(refseqCodingSequence));
            } else {
                refseqCodingSequence = foundRefSeqCodingSequences.get(0);
            }
            logger.info(refseqCodingSequence.toString());

            String location = codingSequenceFeature.getLocation();
            Integer start = null;
            Integer stop = null;

            Matcher m = locationPattern.matcher(location);
            if (m.find()) {
                String startValue = m.group("start");
                start = Integer.valueOf(startValue);
                stop = start;
                String stopValue = m.group("stop");
                if (StringUtils.isNotEmpty(stopValue)) {
                    stop = Integer.valueOf(stopValue);
                }
            }

            if (start != null && stop != null) {

                RegionGroup regionGroup = new RegionGroup();
                regionGroup.setTranscript(transcript);
                regionGroup.setGroupingType(singleGroupingType);
                regionGroup.setId(canvasDAOBeanService.getRegionGroupDAO().save(regionGroup));

                RegionGroupRegionPK rgrPK = new RegionGroupRegionPK(start, stop, "exact", "exact", regionGroup.getId());
                RegionGroupRegion rgr = new RegionGroupRegion(rgrPK);
                logger.info(rgr.toString());
                rgr.setRegionGroup(regionGroup);
                canvasDAOBeanService.getRegionGroupRegionDAO().save(rgr);

                refseqCodingSequence.getLocations().add(regionGroup);
                regionGroup.getRefSeqCodingSequence().add(refseqCodingSequence);
                canvasDAOBeanService.getRefSeqCodingSequenceDAO().save(refseqCodingSequence);

            }

            if (StringUtils.isNotEmpty(ecNumber)) {
                CDSECNumberPK key = new CDSECNumberPK(refseqCodingSequence.getId(), ecNumber);
                CDSECNumber foundCDSECNumber = canvasDAOBeanService.getCDSECNumberDAO().findById(key);
                if (foundCDSECNumber == null) {
                    CDSECNumber cdsECNumber = new CDSECNumber(key);
                    cdsECNumber.setRefseqCodingSequence(refseqCodingSequence);
                    logger.info(cdsECNumber.toString());
                    canvasDAOBeanService.getCDSECNumberDAO().save(cdsECNumber);
                }
            }

            if (CollectionUtils.isNotEmpty(translationExceptionList)) {

                for (TranslationException te : translationExceptionList) {
                    CDSTranslationExceptionPK key = new CDSTranslationExceptionPK(refseqCodingSequence.getId(), te.getRange().getMinimum());
                    CDSTranslationException foundCDSTranslationException = canvasDAOBeanService.getCDSTranslationExceptionDAO()
                            .findById(key);
                    if (foundCDSTranslationException == null) {
                        CDSTranslationException cdsTranslationException = new CDSTranslationException(key);
                        cdsTranslationException.setRefseqCodingSequence(refseqCodingSequence);
                        cdsTranslationException.setStopLocation(te.getRange().getMaximum());
                        cdsTranslationException.setAminoAcid(te.getAminoAcid());
                        logger.info(cdsTranslationException.toString());
                        canvasDAOBeanService.getCDSTranslationExceptionDAO().save(cdsTranslationException);
                    }

                }

            }

        }
    }

    private void persistGenes(String refseqVersion, Transcript transcript, GroupingType singleGroupingType,
            List<org.renci.gbff.model.Feature> features) throws CANVASDAOException {
        logger.debug("ENTERING persistCodingSequence(String, Transcript, GroupingType, List<Feature>)");

        Optional<org.renci.gbff.model.Feature> optionalGeneFeature = features.stream().filter(a -> "gene".equals(a.getType())).findAny();
        if (optionalGeneFeature.isPresent()) {
            org.renci.gbff.model.Feature geneFeature = optionalGeneFeature.get();
            logger.info(geneFeature.toString());

            Map<String, String> qualifiers = geneFeature.getQualifiers();

            String geneName = qualifiers.get("gene");
            String geneDesc = qualifiers.get("note");
            String geneSynonyms = qualifiers.get("gene_synonym");

            RefSeqGene refseqGene = new RefSeqGene(refseqVersion, geneName, geneDesc);

            List<RefSeqGene> foundRefseqGenes = canvasDAOBeanService.getRefSeqGeneDAO().findByExample(refseqGene);
            if (CollectionUtils.isEmpty(foundRefseqGenes)) {
                refseqGene.setId(canvasDAOBeanService.getRefSeqGeneDAO().save(refseqGene));
            } else {
                refseqGene = foundRefseqGenes.get(0);
            }
            logger.info(refseqGene.toString());

            String location = geneFeature.getLocation();
            Integer start = null;
            Integer stop = null;

            Matcher m = locationPattern.matcher(location);
            m.find();
            if (m.matches()) {
                String startValue = m.group("start");
                start = Integer.valueOf(startValue);
                stop = start;
                String stopValue = m.group("stop");
                if (StringUtils.isNotEmpty(stopValue)) {
                    stop = Integer.valueOf(stopValue);
                }
            }

            if (start != null && stop != null) {

                RegionGroup regionGroup = new RegionGroup();
                regionGroup.setTranscript(transcript);
                regionGroup.setGroupingType(singleGroupingType);

                Set<RegionGroup> regionGroups;
                List<RegionGroup> foundRegionGroups = canvasDAOBeanService.getRegionGroupDAO()
                        .findByTranscriptIdAndGroupingType(transcript.getId(), singleGroupingType.getId());
                if (CollectionUtils.isEmpty(foundRegionGroups)) {
                    regionGroup.setId(canvasDAOBeanService.getRegionGroupDAO().save(regionGroup));
                    regionGroups = transcript.getRegionGroups();
                } else {
                    regionGroup = foundRegionGroups.get(0);
                    regionGroups = new HashSet<>(foundRegionGroups);
                }
                logger.info(regionGroup.toString());

                transcript.setRegionGroups(regionGroups);
                canvasDAOBeanService.getTranscriptDAO().save(transcript);

                RegionGroupRegionPK rgrPK = new RegionGroupRegionPK(start, stop, "exact", "exact", regionGroup.getId());

                RegionGroupRegion foundRGR = canvasDAOBeanService.getRegionGroupRegionDAO().findById(rgrPK);
                if (foundRGR == null) {
                    RegionGroupRegion rgr = new RegionGroupRegion(rgrPK);
                    rgr.setRegionGroup(regionGroup);
                    logger.info(rgr.toString());
                    canvasDAOBeanService.getRegionGroupRegionDAO().save(rgr);
                }

                refseqGene.getLocations().add(regionGroup);
                regionGroup.getRefSeqGenes().add(refseqGene);
                canvasDAOBeanService.getRefSeqGeneDAO().save(refseqGene);

            }

            AnnotationGene annotationGene = null;

            List<AnnotationGene> foundAnnotationGenes = canvasDAOBeanService.getAnnotationGeneDAO().findByName(geneName);
            if (CollectionUtils.isEmpty(foundAnnotationGenes)) {
                annotationGene = new AnnotationGene(geneName, geneDesc);
                annotationGene.setId(canvasDAOBeanService.getAnnotationGeneDAO().save(annotationGene));
            } else {
                annotationGene = foundAnnotationGenes.get(0);
            }
            logger.info(annotationGene.toString());

            Set<AnnotationGeneExternalId> externalIdSet = new HashSet<>();

            AnnotationGeneExternalIdPK annotationGeneExternalIdPK = new AnnotationGeneExternalIdPK(refseqGene.getId(),
                    annotationGene.getId(), "refseq", refseqVersion);

            AnnotationGeneExternalId foundAnnotationGeneExternalId = canvasDAOBeanService.getAnnotationGeneExternalIdDAO()
                    .findById(annotationGeneExternalIdPK);

            if (foundAnnotationGeneExternalId == null) {

                AnnotationGeneExternalId annotationGeneExternalId = new AnnotationGeneExternalId(annotationGeneExternalIdPK);

                annotationGeneExternalId.setGene(annotationGene);
                logger.info(annotationGeneExternalId.toString());
                canvasDAOBeanService.getAnnotationGeneExternalIdDAO().save(annotationGeneExternalId);
                externalIdSet.add(annotationGeneExternalId);

            } else {
                externalIdSet.add(foundAnnotationGeneExternalId);
            }

            annotationGene.setExternals(externalIdSet);
            canvasDAOBeanService.getAnnotationGeneDAO().save(annotationGene);

            List<AnnotationGeneSynonym> annotationGeneSynonyms = canvasDAOBeanService.getAnnotationGeneSynonymDAO()
                    .findByGeneId(annotationGene.getId());

            Set<AnnotationGeneSynonym> annotationGeneSynonymSet = new HashSet<>(annotationGeneSynonyms);

            List<String> synonymList = new ArrayList<>();
            if (StringUtils.isNotEmpty(geneSynonyms)) {
                if (geneSynonyms.contains(";")) {
                    synonymList.addAll(Arrays.asList(geneSynonyms.split(";")));
                } else {
                    synonymList.add(geneSynonyms);
                }
            }

            if (CollectionUtils.isNotEmpty(synonymList)) {
                for (String synonym : synonymList) {
                    AnnotationGeneSynonymPK annotationGeneSynonymPK = new AnnotationGeneSynonymPK(annotationGene.getId(), synonym);
                    AnnotationGeneSynonym foundAnnotationGeneSynonym = canvasDAOBeanService.getAnnotationGeneSynonymDAO()
                            .findById(annotationGeneSynonymPK);
                    if (foundAnnotationGeneSynonym == null) {
                        AnnotationGeneSynonym annotationGeneSynonym = new AnnotationGeneSynonym(annotationGeneSynonymPK);
                        annotationGeneSynonym.setGene(annotationGene);
                        annotationGeneSynonym.setId(canvasDAOBeanService.getAnnotationGeneSynonymDAO().save(annotationGeneSynonym));
                        logger.info(annotationGeneSynonym.toString());
                        annotationGeneSynonymSet.add(annotationGeneSynonym);
                    } else {
                        logger.info(foundAnnotationGeneSynonym.toString());
                        annotationGeneSynonymSet.add(foundAnnotationGeneSynonym);
                    }
                }
            }

            annotationGene.setSynonyms(annotationGeneSynonymSet);
            canvasDAOBeanService.getAnnotationGeneDAO().save(annotationGene);

        }

    }

    private Transcript persistTranscript(String refseqVersion, Sequence sequence) throws CANVASDAOException {
        logger.debug("ENTERING persistTranscript(String, Sequence)");

        StringBuilder seq = new StringBuilder();
        sequence.getOrigin().forEach(a -> seq.append(a.getSequence()));
        String accession = sequence.getAccession();
        if (sequence.getAccession().contains(" ")) {
            accession = sequence.getAccession().substring(0, sequence.getAccession().indexOf(" "));
        }

        Transcript transcript = new Transcript(sequence.getVersion().trim(), accession, seq.toString());

        Transcript foundTranscript = canvasDAOBeanService.getTranscriptDAO().findById(sequence.getVersion().trim());
        if (foundTranscript == null) {
            canvasDAOBeanService.getTranscriptDAO().save(transcript);
        }

        TranscriptRefSeqVersionPK transcriptRefSeqVersionPK = new TranscriptRefSeqVersionPK(sequence.getVersion().trim(), refseqVersion);
        TranscriptRefSeqVersion foundTranscriptRefSeqVersion = canvasDAOBeanService.getTranscriptRefSeqVersionDAO()
                .findById(transcriptRefSeqVersionPK);
        if (foundTranscriptRefSeqVersion == null) {
            TranscriptRefSeqVersion transcriptRefSeqVersion = new TranscriptRefSeqVersion(transcriptRefSeqVersionPK);
            transcriptRefSeqVersion.setTranscript(transcript);
            logger.info(transcriptRefSeqVersion.toString());
            canvasDAOBeanService.getTranscriptRefSeqVersionDAO().save(transcriptRefSeqVersion);

            transcript.getRefseqVersions().add(transcriptRefSeqVersion);
            canvasDAOBeanService.getTranscriptDAO().save(transcript);
        }

        return transcript;
    }

}
