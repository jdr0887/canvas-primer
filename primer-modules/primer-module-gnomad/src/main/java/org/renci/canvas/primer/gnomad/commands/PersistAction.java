package org.renci.canvas.primer.gnomad.commands;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.Option;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.commons.LocatedVariantFactory;
import org.renci.canvas.dao.gnomad.model.GnomADVariantFrequency;
import org.renci.canvas.dao.gnomad.model.GnomADVariantFrequencyPK;
import org.renci.canvas.dao.ref.model.GenomeRef;
import org.renci.canvas.dao.ref.model.GenomeRefSeq;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.renci.canvas.dao.var.model.VariantType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

@Command(scope = "gnomad", name = "persist", description = "")
@Service
public class PersistAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(PersistAction.class);

    @Reference
    private CANVASDAOBeanService canvasDAOBeanService;

    @Option(name = "--gnomadExomesVCF", required = true, multiValued = false)
    private String gnomadExomesVCF;

    public PersistAction() {
        super();
    }

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");

        Pair<String, File> versionAndFilePair = downloadLatest();

        String version = versionAndFilePair.getLeft();
        File latestFile = versionAndFilePair.getRight();

        Executors.newSingleThreadExecutor().submit(() -> {

            try {
                List<GenomeRef> allGenomeRefs = canvasDAOBeanService.getGenomeRefDAO().findAll();
                GenomeRef genomeRef = allGenomeRefs.stream().filter(a -> a.getName().startsWith("38"))
                        .sorted((a, b) -> a.getName().compareTo(b.getName())).findFirst().get();

                List<VariantType> allVariantTypes = canvasDAOBeanService.getVariantTypeDAO().findAll();

                List<GenomeRefSeq> allGenomeRefSeqs = canvasDAOBeanService.getGenomeRefSeqDAO()
                        .findByGenomeRefIdAndSeqType(genomeRef.getId(), "Chromosome");

                try (VCFFileReader vcfFileReader = new VCFFileReader(latestFile, false)) {

                    for (VariantContext variantContext : vcfFileReader) {

                        GenomeRefSeq genomeRefSeq = allGenomeRefSeqs.stream().filter(a -> a.getContig().equals(variantContext.getContig()))
                                .findFirst().get();

                        CommonInfo commonInfo = variantContext.getCommonInfo();

                        for (Allele altAllele : variantContext.getAlternateAlleles()) {

                            LocatedVariant locatedVariant = LocatedVariantFactory.create(genomeRef, genomeRefSeq, variantContext, altAllele,
                                    allVariantTypes);

                            logger.info(locatedVariant.toString());
                            List<LocatedVariant> foundLocatedVariants = canvasDAOBeanService.getLocatedVariantDAO()
                                    .findByExample(locatedVariant);
                            if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                                locatedVariant = foundLocatedVariants.get(0);
                            } else {
                                locatedVariant.setId(canvasDAOBeanService.getLocatedVariantDAO().save(locatedVariant));
                            }

                            final LocatedVariant finalLocatedVariant = locatedVariant;

                            int alleleIndex = variantContext.getAlleleIndex(altAllele);
                            logger.debug("alleleIndex: {}", alleleIndex);

                            Arrays.asList("AFR", "AMR", "ASJ", "EAS", "FIN", "NFE", "OTH", "SAS").parallelStream().forEach(population -> {

                                try {

                                    String alleleFrequencyValue = commonInfo.getAttributeAsString(String.format("AF_%s", population), "");
                                    String alleleCountValue = commonInfo.getAttributeAsString(String.format("AC_%s", population), "");
                                    String alleleTotalValue = commonInfo.getAttributeAsString(String.format("AN_%s", population), "");
                                    String hemizygousCountValue = commonInfo.getAttributeAsString(String.format("Hemi_%s", population), "");
                                    String homozygousCountValue = commonInfo.getAttributeAsString(String.format("Hom_%s", population), "");

                                    GnomADVariantFrequencyPK variantFrequencyPK = new GnomADVariantFrequencyPK(finalLocatedVariant.getId(),
                                            version, population);

                                    GnomADVariantFrequency foundGnomADVariantFrequency = canvasDAOBeanService.getGnomADVariantFrequencyDAO()
                                            .findById(variantFrequencyPK);

                                    GnomADVariantFrequency variantFrequency = null;
                                    if (foundGnomADVariantFrequency == null) {
                                        variantFrequency = new GnomADVariantFrequency(variantFrequencyPK);
                                    } else {
                                        variantFrequency = foundGnomADVariantFrequency;
                                    }
                                    variantFrequency.setLocatedVariant(finalLocatedVariant);

                                    Integer alleleTotal = StringUtils.isNotEmpty(alleleTotalValue) && !".".equals(alleleTotalValue)
                                            ? Integer.valueOf(alleleTotalValue) : 0;
                                    Double alleleFrequency = null;
                                    Integer alleleCount = null;
                                    Integer hemizygousCount = null;
                                    Integer homozygousCount = null;

                                    if (alleleIndex > 0) {

                                        List<String> alleleFrequencyValues = commonInfo
                                                .getAttributeAsStringList(String.format("AF_%s", population), "");
                                        if (CollectionUtils.isNotEmpty(alleleFrequencyValues)) {
                                            alleleFrequencyValue = alleleFrequencyValues.get(alleleIndex - 1);
                                            alleleFrequency = StringUtils.isNotEmpty(alleleFrequencyValue)
                                                    && !".".equals(alleleFrequencyValue) ? Double.valueOf(alleleFrequencyValue) : 0D;
                                        } else {
                                            alleleFrequency = 0D;
                                        }

                                        List<String> alleleCountValues = commonInfo
                                                .getAttributeAsStringList(String.format("AC_%s", population), "");
                                        if (CollectionUtils.isNotEmpty(alleleCountValues)) {
                                            alleleCountValue = alleleCountValues.get(alleleIndex - 1);
                                            alleleCount = StringUtils.isNotEmpty(alleleCountValue) && !".".equals(alleleCountValue)
                                                    ? Integer.valueOf(alleleCountValue) : 0;
                                        } else {
                                            alleleCount = 0;
                                        }

                                        List<String> hemizygousCountValues = commonInfo
                                                .getAttributeAsStringList(String.format("Hemi_%s", population), "");
                                        if (CollectionUtils.isNotEmpty(hemizygousCountValues)) {
                                            hemizygousCountValue = hemizygousCountValues.get(alleleIndex - 1);
                                            hemizygousCount = StringUtils.isNotEmpty(hemizygousCountValue)
                                                    && !".".equals(hemizygousCountValue) ? Integer.valueOf(hemizygousCountValue) : 0;
                                        } else {
                                            hemizygousCount = 0;
                                        }

                                        List<String> homozygousCountValues = commonInfo
                                                .getAttributeAsStringList(String.format("Hom_%s", population), "");
                                        if (CollectionUtils.isNotEmpty(homozygousCountValues)) {
                                            homozygousCountValue = homozygousCountValues.get(alleleIndex - 1);
                                            homozygousCount = StringUtils.isNotEmpty(homozygousCountValue)
                                                    && !".".equals(homozygousCountValue) ? Integer.valueOf(homozygousCountValue) : 0;
                                        } else {
                                            homozygousCount = 0;
                                        }

                                    } else {

                                        alleleFrequency = StringUtils.isNotEmpty(alleleFrequencyValue) && !".".equals(alleleFrequencyValue)
                                                ? Double.valueOf(alleleFrequencyValue) : 0D;
                                        alleleCount = StringUtils.isNotEmpty(alleleCountValue) && !".".equals(alleleCountValue)
                                                ? Integer.valueOf(alleleCountValue) : 0;
                                        alleleTotal = StringUtils.isNotEmpty(alleleTotalValue) && !".".equals(alleleTotalValue)
                                                ? Integer.valueOf(alleleTotalValue) : 0;
                                        hemizygousCount = StringUtils.isNotEmpty(hemizygousCountValue) && !".".equals(hemizygousCountValue)
                                                ? Integer.valueOf(hemizygousCountValue) : 0;
                                        homozygousCount = StringUtils.isNotEmpty(homozygousCountValue) && !".".equals(homozygousCountValue)
                                                ? Integer.valueOf(homozygousCountValue) : 0;

                                    }

                                    variantFrequency.setAlternateAlleleFrequency(alleleFrequency);
                                    variantFrequency.setAlternateAlleleCount(alleleCount);
                                    variantFrequency.setTotalAlleleCount(alleleTotal);
                                    variantFrequency.setHemizygousCount(hemizygousCount);
                                    variantFrequency.setHomozygousCount(homozygousCount);

                                    canvasDAOBeanService.getGnomADVariantFrequencyDAO().save(variantFrequency);

                                    logger.debug(variantFrequency.toString());
                                } catch (NumberFormatException | CANVASDAOException e) {
                                    logger.error(e.getMessage(), e);
                                }

                            });

                        }

                    }
                }
            } catch (Exception e) {
                logger.error(e.getMessage(), e);
            }

        });

        return null;
    }

    private Pair<String, File> downloadLatest() {

        File download = new File(getGnomadExomesVCF());

        Pattern p = Pattern.compile("gnomad\\.exomes\\.r(?<version>\\d\\.\\d\\.\\d)\\.sites\\.vcf");
        Matcher m = p.matcher(download.getName());
        m.find();
        String version = null;
        if (m.matches()) {
            version = m.group("version");
        }

        return Pair.of(version, download);
    }

    public String getGnomadExomesVCF() {
        return gnomadExomesVCF;
    }

    public void setGnomadExomesVCF(String gnomadExomesVCF) {
        this.gnomadExomesVCF = gnomadExomesVCF;
    }

}
