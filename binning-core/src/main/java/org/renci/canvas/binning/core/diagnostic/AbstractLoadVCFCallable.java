package org.renci.canvas.binning.core.diagnostic;

import java.io.File;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.ResourceBundle;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.renci.canvas.binning.core.BinningException;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.DiagnosticBinningJob;
import org.renci.canvas.dao.ref.model.GenomeRef;
import org.renci.canvas.dao.ref.model.GenomeRefSeq;
import org.renci.canvas.dao.var.model.Assembly;
import org.renci.canvas.dao.var.model.AssemblyLocatedVariant;
import org.renci.canvas.dao.var.model.AssemblyLocatedVariantPK;
import org.renci.canvas.dao.var.model.AssemblyLocatedVariantQC;
import org.renci.canvas.dao.var.model.AssemblyLocatedVariantQCPK;
import org.renci.canvas.dao.var.model.CanonicalAllele;
import org.renci.canvas.dao.var.model.Lab;
import org.renci.canvas.dao.var.model.Library;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.renci.canvas.dao.var.model.Project;
import org.renci.canvas.dao.var.model.Sample;
import org.renci.canvas.dao.var.model.VariantSet;
import org.renci.canvas.dao.var.model.VariantSetLoad;
import org.renci.canvas.dao.var.model.VariantType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public abstract class AbstractLoadVCFCallable implements Callable<Void> {

    private static final Logger logger = LoggerFactory.getLogger(AbstractLoadVCFCallable.class);

    private CANVASDAOBeanService daoBean;

    private DiagnosticBinningJob binningJob;

    public abstract File getVCF(String participant) throws BinningException;

    public abstract Set<String> getExcludesFilter();

    public abstract String getLabName();

    public abstract String getLibraryName();

    public abstract String getStudyName();

    public abstract GenomeRef getGenomeRef();

    public abstract LocatedVariant liftOver(LocatedVariant locatedVariant) throws BinningException;

    public AbstractLoadVCFCallable(CANVASDAOBeanService daoBean, DiagnosticBinningJob binningJob) {
        super();
        this.daoBean = daoBean;
        this.binningJob = binningJob;
    }

    @Override
    public Void call() throws BinningException {
        logger.debug("ENTERING run()");

        int refSkippedCount = 0;
        int errorCount = 0;
        int filteredCount = 0;
        int noCallCount = 0;

        try {

            Date startDate = new Date();

            List<VariantType> allVariantTypes = daoBean.getVariantTypeDAO().findAll();

            File vcfFile = null;
            if (StringUtils.isEmpty(binningJob.getVcfFile())) {
                vcfFile = getVCF(binningJob.getParticipant());
                binningJob.setVcfFile(vcfFile.getAbsolutePath());
                daoBean.getDiagnosticBinningJobDAO().save(binningJob);
            } else {
                vcfFile = new File(binningJob.getVcfFile());
                if (!vcfFile.exists()) {
                    vcfFile = getVCF(binningJob.getParticipant());
                    binningJob.setVcfFile(vcfFile.getAbsolutePath());
                    daoBean.getDiagnosticBinningJobDAO().save(binningJob);
                }
            }
            logger.info(binningJob.toString());

            final GenomeRef genomeRef = getGenomeRef();
            logger.info(genomeRef.toString());

            Map<String, List<VariantContext>> variantContext2SampleNameMap = new HashMap<String, List<VariantContext>>();

            Set<String> excludesFilter = getExcludesFilter();
            int expectecedLocatedVariantCount = 0;

            try (VCFFileReader vcfFileReader = new VCFFileReader(vcfFile, false)) {

                VCFHeader vcfHeader = vcfFileReader.getFileHeader();
                List<String> sampleNames = vcfHeader.getGenotypeSamples();

                for (String sampleName : sampleNames) {

                    if (!variantContext2SampleNameMap.containsKey(sampleName)) {
                        variantContext2SampleNameMap.put(sampleName, new ArrayList<VariantContext>());
                    }

                    variantContextLoop: for (VariantContext variantContext : vcfFileReader) {

                        Allele refAllele = variantContext.getReference();

                        if (refAllele.isNoCall()) {
                            noCallCount++;
                            continue;
                        }

                        String altAlleles = StringUtils.join(variantContext.getAlternateAlleles().toArray());
                        if (!altAlleles.matches("[AaCcGgTt,]*")) {
                            errorCount++;
                            continue;
                        }

                        if (CollectionUtils.containsAny(variantContext.getFilters(), excludesFilter)) {
                            filteredCount++;
                            continue;
                        }

                        GenotypesContext genotypesContext = variantContext.getGenotypes();

                        for (Genotype genotype : genotypesContext) {

                            if (genotype.isNoCall()) {
                                noCallCount++;
                                continue variantContextLoop;
                            }

                            if (genotype.isHomRef()) {
                                refSkippedCount++;
                                continue variantContextLoop;
                            }
                        }

                        for (Allele altAllele : variantContext.getAlternateAlleles()) {
                            expectecedLocatedVariantCount++;
                        }

                        variantContext2SampleNameMap.get(sampleName).add(variantContext);

                    }

                }
            } catch (Exception e) {
                logger.error(e.getMessage(), e);
                throw new BinningException(e);
            }
            logger.info("expectecedLocatedVariantCount: {}", expectecedLocatedVariantCount);

            final List<GenomeRefSeq> allGenomeRefSeqs = daoBean.getGenomeRefSeqDAO().findByGenomeRefIdAndSeqType(genomeRef.getId(),
                    "Chromosome");

            final Set<LocatedVariant> locatedVariantSet = new HashSet<>();

            for (String sampleName : variantContext2SampleNameMap.keySet()) {

                if (!sampleName.equals(binningJob.getParticipant())) {
                    throw new BinningException("sampleName does match participant");
                }
                Assembly assembly = buildAssembly(sampleName);
                logger.info(assembly.toString());

                List<VariantContext> variantContextList = variantContext2SampleNameMap.get(sampleName);

                if (CollectionUtils.isNotEmpty(variantContextList)) {

                    ExecutorService es = Executors.newFixedThreadPool(6);

                    // snp
                    for (VariantContext variantContext : variantContextList) {

                        for (Allele altAllele : variantContext.getAlternateAlleles()) {

                            if (variantContext.isSNP()) {

                                es.submit(() -> {

                                    Optional<GenomeRefSeq> genomeRefSeqOptional = allGenomeRefSeqs.stream()
                                            .filter(a -> a.getId().equals(variantContext.getContig())).findAny();
                                    if (genomeRefSeqOptional.isPresent()) {

                                        GenomeRefSeq genomeRefSeq = genomeRefSeqOptional.get();

                                        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);
                                        locatedVariant.setSeq(altAllele.getDisplayString());
                                        locatedVariant.setRef(variantContext.getReference().getDisplayString());
                                        locatedVariant.setPosition(variantContext.getStart());
                                        locatedVariant.setVariantType(
                                                allVariantTypes.stream().filter(a -> a.getId().equals("snp")).findAny().get());
                                        locatedVariant.setEndPosition(variantContext.getStart() + locatedVariant.getRef().length());

                                        try {
                                            List<LocatedVariant> foundLocatedVariants = daoBean.getLocatedVariantDAO()
                                                    .findByExample(locatedVariant);
                                            if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                                                locatedVariant = foundLocatedVariants.get(0);
                                            } else {
                                                locatedVariant.setId(daoBean.getLocatedVariantDAO().save(locatedVariant));
                                            }
                                            logger.info(locatedVariant.toString());
                                            createAssmeblyLocatedVariantQC(sampleName, variantContext, locatedVariant, assembly);
                                            locatedVariantSet.add(locatedVariant);
                                        } catch (CANVASDAOException | BinningException e) {
                                            logger.error(e.getMessage(), e);
                                        }
                                    }

                                });

                            }
                        }

                    }
                    es.shutdown();
                    if (!es.awaitTermination(1L, TimeUnit.HOURS)) {
                        es.shutdownNow();
                    }

                    es = Executors.newFixedThreadPool(6);
                    // ins
                    for (VariantContext variantContext : variantContextList) {

                        for (Allele altAllele : variantContext.getAlternateAlleles()) {

                            if (variantContext.isIndel() && variantContext.isSimpleInsertion()) {

                                es.submit(() -> {

                                    Optional<GenomeRefSeq> genomeRefSeqOptional = allGenomeRefSeqs.stream()
                                            .filter(a -> a.getId().equals(variantContext.getContig())).findAny();
                                    if (genomeRefSeqOptional.isPresent()) {
                                        GenomeRefSeq genomeRefSeq = genomeRefSeqOptional.get();

                                        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);
                                        locatedVariant.setPosition(variantContext.getStart());
                                        locatedVariant.setVariantType(
                                                allVariantTypes.stream().filter(a -> a.getId().equals("ins")).findAny().get());
                                        String ref = variantContext.getReference().getDisplayString();
                                        locatedVariant.setSeq(altAllele.getDisplayString().replaceFirst(ref, ""));
                                        locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
                                        locatedVariant.setRef("");

                                        try {
                                            List<LocatedVariant> foundLocatedVariants = daoBean.getLocatedVariantDAO()
                                                    .findByExample(locatedVariant);
                                            if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                                                locatedVariant = foundLocatedVariants.get(0);
                                            } else {
                                                locatedVariant.setId(daoBean.getLocatedVariantDAO().save(locatedVariant));
                                            }
                                            logger.info(locatedVariant.toString());
                                            createAssmeblyLocatedVariantQC(sampleName, variantContext, locatedVariant, assembly);
                                            locatedVariantSet.add(locatedVariant);
                                        } catch (CANVASDAOException | BinningException e) {
                                            logger.error(e.getMessage(), e);
                                        }
                                    }

                                });

                            }
                        }

                    }
                    es.shutdown();
                    if (!es.awaitTermination(1L, TimeUnit.HOURS)) {
                        es.shutdownNow();
                    }

                    es = Executors.newFixedThreadPool(6);
                    // del
                    for (VariantContext variantContext : variantContextList) {

                        for (Allele altAllele : variantContext.getAlternateAlleles()) {

                            if (variantContext.isIndel() && variantContext.isSimpleDeletion()) {

                                es.submit(() -> {

                                    Optional<GenomeRefSeq> genomeRefSeqOptional = allGenomeRefSeqs.parallelStream()
                                            .filter(a -> a.getId().equals(variantContext.getContig())).findAny();
                                    if (genomeRefSeqOptional.isPresent()) {

                                        GenomeRefSeq genomeRefSeq = genomeRefSeqOptional.get();
                                        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);

                                        locatedVariant.setPosition(variantContext.getStart() + 1);
                                        locatedVariant.setRef(variantContext.getReference().getDisplayString()
                                                .replaceFirst(altAllele.getDisplayString(), ""));
                                        locatedVariant.setSeq(locatedVariant.getRef());
                                        locatedVariant.setVariantType(
                                                allVariantTypes.stream().filter(a -> a.getId().equals("del")).findAny().get());
                                        locatedVariant.setEndPosition(locatedVariant.getPosition() + locatedVariant.getRef().length());

                                        try {
                                            List<LocatedVariant> foundLocatedVariants = daoBean.getLocatedVariantDAO()
                                                    .findByExample(locatedVariant);
                                            if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                                                locatedVariant = foundLocatedVariants.get(0);
                                            } else {
                                                locatedVariant.setId(daoBean.getLocatedVariantDAO().save(locatedVariant));
                                            }
                                            logger.info(locatedVariant.toString());
                                            createAssmeblyLocatedVariantQC(sampleName, variantContext, locatedVariant, assembly);
                                            locatedVariantSet.add(locatedVariant);
                                        } catch (CANVASDAOException | BinningException e) {
                                            logger.error(e.getMessage(), e);
                                        }

                                    }

                                });

                            }
                        }

                    }
                    es.shutdown();
                    if (!es.awaitTermination(1L, TimeUnit.HOURS)) {
                        es.shutdownNow();
                    }

                    es = Executors.newFixedThreadPool(6);

                    for (VariantContext variantContext : variantContextList) {

                        Optional<GenomeRefSeq> genomeRefSeqOptional = allGenomeRefSeqs.parallelStream()
                                .filter(a -> a.getId().equals(variantContext.getContig())).findAny();
                        if (genomeRefSeqOptional.isPresent()) {
                            GenomeRefSeq genomeRefSeq = genomeRefSeqOptional.get();

                            // cant trust htsjdk to parse properly...switch on freebayes type (if available)
                            List<String> types = variantContext.getAttributeAsStringList("TYPE", null);
                            for (Allele altAllele : variantContext.getAlternateAlleles()) {

                                if (CollectionUtils.isNotEmpty(types) && !variantContext.isSNP()
                                        && !(variantContext.isIndel() && variantContext.isSimpleInsertion())
                                        && !(variantContext.isIndel() && variantContext.isSimpleDeletion())
                                        && (variantContext.isComplexIndel() || variantContext.isMNP())) {

                                    es.submit(() -> {
                                        String ref = variantContext.getReference().getDisplayString();
                                        String alt = altAllele.getDisplayString();

                                        char[] referenceChars = ref.toCharArray();
                                        char[] alternateChars = alt.toCharArray();

                                        LocatedVariant locatedVariant = new LocatedVariant(genomeRef, genomeRefSeq);
                                        String type = types.get(variantContext.getAlleleIndex(altAllele) - 1);
                                        switch (type) {
                                            case "del":
                                                locatedVariant.setVariantType(
                                                        allVariantTypes.stream().filter(a -> a.getId().equals("del")).findAny().get());
                                                locatedVariant.setPosition(variantContext.getStart() + 1);
                                                locatedVariant.setEndPosition(variantContext.getStart() + 1 + ref.length());
                                                locatedVariant.setRef(ref);
                                                locatedVariant.setSeq(ref);
                                                break;
                                            case "ins":
                                                locatedVariant.setVariantType(
                                                        allVariantTypes.stream().filter(a -> a.getId().equals("ins")).findAny().get());
                                                locatedVariant.setPosition(variantContext.getStart());

                                                if (referenceChars.length > 1 && alternateChars.length > 1) {

                                                    StringBuilder frontChars2Remove = new StringBuilder();
                                                    StringBuilder backChars2Remove = new StringBuilder();

                                                    for (int i = 0; i < referenceChars.length; ++i) {
                                                        if (referenceChars[i] != alternateChars[i]) {
                                                            break;
                                                        }
                                                        frontChars2Remove.append(referenceChars[i]);
                                                    }

                                                    for (int i = referenceChars.length - 1; i > 0; --i) {
                                                        if (referenceChars[i] != alternateChars[i]) {
                                                            break;
                                                        }
                                                        backChars2Remove.append(referenceChars[i]);
                                                    }

                                                    if (frontChars2Remove.length() > 0) {
                                                        ref = ref.replaceFirst(frontChars2Remove.toString(), "");
                                                        alt = alt.replaceFirst(frontChars2Remove.toString(), "");
                                                    }

                                                    if (backChars2Remove.length() > 0) {
                                                        backChars2Remove.reverse();
                                                        ref = StringUtils.removeEnd(ref, backChars2Remove.toString());
                                                        alt = StringUtils.removeEnd(alt, backChars2Remove.toString());
                                                    }

                                                    locatedVariant.setPosition(variantContext.getStart()
                                                            + (frontChars2Remove.length() > 0 ? frontChars2Remove.length() : 0));
                                                }

                                                locatedVariant.setEndPosition(locatedVariant.getPosition() + 1);
                                                locatedVariant.setRef("");
                                                locatedVariant.setSeq(alt);

                                                break;
                                            case "snp":
                                                locatedVariant.setVariantType(
                                                        allVariantTypes.stream().filter(a -> a.getId().equals("snp")).findAny().get());
                                                locatedVariant.setPosition(variantContext.getStart());

                                                if (referenceChars.length > 1 && alternateChars.length > 1) {

                                                    StringBuilder frontChars2Remove = new StringBuilder();
                                                    StringBuilder backChars2Remove = new StringBuilder();

                                                    for (int i = 0; i < referenceChars.length; ++i) {
                                                        if (referenceChars[i] != alternateChars[i]) {
                                                            break;
                                                        }
                                                        frontChars2Remove.append(referenceChars[i]);
                                                    }

                                                    for (int i = referenceChars.length - 1; i > 0; --i) {
                                                        if (referenceChars[i] != alternateChars[i]) {
                                                            break;
                                                        }
                                                        backChars2Remove.append(referenceChars[i]);
                                                    }

                                                    if (frontChars2Remove.length() > 0) {
                                                        ref = ref.replaceFirst(frontChars2Remove.toString(), "");
                                                        alt = alt.replaceFirst(frontChars2Remove.toString(), "");
                                                    }

                                                    if (backChars2Remove.length() > 0) {
                                                        backChars2Remove.reverse();
                                                        ref = StringUtils.removeEnd(ref, backChars2Remove.toString());
                                                        alt = StringUtils.removeEnd(alt, backChars2Remove.toString());
                                                    }

                                                    locatedVariant.setPosition(variantContext.getStart()
                                                            + (frontChars2Remove.length() > 0 ? frontChars2Remove.length() : 0));
                                                }
                                                locatedVariant.setEndPosition(locatedVariant.getPosition() + 1);
                                                locatedVariant.setRef(ref);
                                                locatedVariant.setSeq(alt);

                                                break;
                                            default:
                                                locatedVariant.setVariantType(
                                                        allVariantTypes.stream().filter(a -> a.getId().equals("sub")).findAny().get());
                                                locatedVariant.setPosition(variantContext.getStart());
                                                locatedVariant.setRef(ref);
                                                locatedVariant.setSeq(alt);
                                                locatedVariant.setEndPosition(locatedVariant.getPosition() + ref.length());
                                                break;
                                        }

                                        try {
                                            List<LocatedVariant> foundLocatedVariants = daoBean.getLocatedVariantDAO()
                                                    .findByExample(locatedVariant);
                                            if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                                                locatedVariant = foundLocatedVariants.get(0);
                                            } else {
                                                locatedVariant.setId(daoBean.getLocatedVariantDAO().save(locatedVariant));
                                            }

                                            if (!locatedVariantSet.contains(locatedVariant)) {
                                                logger.info(locatedVariant.toString());
                                                createAssmeblyLocatedVariantQC(sampleName, variantContext, locatedVariant, assembly);
                                                locatedVariantSet.add(locatedVariant);
                                            }
                                        } catch (CANVASDAOException | BinningException e) {
                                            logger.error(e.getMessage(), e);
                                        }

                                    });

                                }

                            }

                        }

                    }
                    es.shutdown();
                    if (!es.awaitTermination(1L, TimeUnit.HOURS)) {
                        es.shutdownNow();
                    }
                }

                ExecutorService es = Executors.newFixedThreadPool(4);
                for (LocatedVariant locatedVariant : locatedVariantSet) {
                    es.submit(() -> {
                        try {
                            logger.info(locatedVariant.toString());
                            canonicalize(locatedVariant);
                        } catch (Exception e) {
                            logger.error(e.getMessage(), e);
                        }
                    });
                }
                es.shutdown();
                if (!es.awaitTermination(2L, TimeUnit.HOURS)) {
                    es.shutdownNow();
                }

                VariantSetLoad variantSetLoad = new VariantSetLoad();
                variantSetLoad.setLoadFilename(vcfFile.getAbsolutePath());
                variantSetLoad.setLoadProgramName(getClass().getName());

                // BundleContext bundleContext = FrameworkUtil.getBundle(getClass()).getBundleContext();
                // Bundle bundle = bundleContext.getBundle();
                // String version = bundle.getVersion().toString();
                String version = ResourceBundle.getBundle("org/renci/canvas/binning/binning").getString("version");

                variantSetLoad.setLoadProgramVersion(version);
                variantSetLoad.setId(assembly.getVariantSet().getId());
                variantSetLoad.setVariantSet(assembly.getVariantSet());

                List<VariantSetLoad> foundVariantSetLoads = daoBean.getVariantSetLoadDAO().findByExample(variantSetLoad);
                if (CollectionUtils.isNotEmpty(foundVariantSetLoads)) {
                    variantSetLoad = foundVariantSetLoads.get(0);
                }

                variantSetLoad.setLoadUser(System.getProperty("user.name"));
                variantSetLoad.setLoadTimeStart(startDate);
                variantSetLoad.setLoadTimeStop(new Date());
                variantSetLoad.setNotes("");
                variantSetLoad.setNumberOfDelRows(locatedVariantSet.stream().filter(a -> a.getVariantType().getId().equals("del"))
                        .collect(Collectors.toList()).size());
                variantSetLoad.setNumberOfInsRows(locatedVariantSet.stream().filter(a -> a.getVariantType().getId().equals("ins"))
                        .collect(Collectors.toList()).size());
                variantSetLoad.setNumberOfSNPRows(locatedVariantSet.stream().filter(a -> a.getVariantType().getId().equals("snp"))
                        .collect(Collectors.toList()).size());
                variantSetLoad.setNumberOfSubRows(locatedVariantSet.stream().filter(a -> a.getVariantType().getId().equals("sub"))
                        .collect(Collectors.toList()).size());
                variantSetLoad.setNumberOfErrorRows(errorCount);
                variantSetLoad.setNumberOfFilteredRows(filteredCount);
                variantSetLoad.setNumberOfMultiRows(0);
                variantSetLoad.setNumberOfSkippedRefRows(refSkippedCount);
                daoBean.getVariantSetLoadDAO().save(variantSetLoad);
                logger.info(variantSetLoad.toString());

            }

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            throw new BinningException(e);
        }
        return null;
    }

    private void canonicalize(LocatedVariant locatedVariant) throws BinningException {
        logger.debug("ENTERING canonicalize(LocatedVariant)");
        try {

            CanonicalAllele canonicalAllele = null;
            // first try to find CanonicalAllele by LocatedVariant
            List<CanonicalAllele> foundCanonicalAlleles = daoBean.getCanonicalAlleleDAO().findByLocatedVariantId(locatedVariant.getId());
            if (CollectionUtils.isNotEmpty(foundCanonicalAlleles)) {
                canonicalAllele = foundCanonicalAlleles.get(0);
            }

            LocatedVariant liftOverLocatedVariant = liftOver(locatedVariant);
            if (liftOverLocatedVariant != null) {
                if (locatedVariant.getVariantType().getId().equals("ins")) {
                    // could have had a deletion in ref
                    liftOverLocatedVariant.setEndPosition(liftOverLocatedVariant.getPosition() + 1);
                }
                List<LocatedVariant> foundLocatedVariants = daoBean.getLocatedVariantDAO().findByExample(liftOverLocatedVariant);
                if (CollectionUtils.isNotEmpty(foundLocatedVariants)) {
                    liftOverLocatedVariant = foundLocatedVariants.get(0);
                } else {
                    liftOverLocatedVariant.setId(daoBean.getLocatedVariantDAO().save(liftOverLocatedVariant));
                }
                logger.info("liftOver: {}", liftOverLocatedVariant.toString());
            }

            // if not found, try to find CanonicalAllele by liftover LocatedVariant
            if (liftOverLocatedVariant != null && canonicalAllele == null) {
                List<CanonicalAllele> foundCanonicalAllelesByLiftOverLocatedVariant = daoBean.getCanonicalAlleleDAO()
                        .findByLocatedVariantId(liftOverLocatedVariant.getId());
                if (CollectionUtils.isNotEmpty(foundCanonicalAllelesByLiftOverLocatedVariant)) {
                    canonicalAllele = foundCanonicalAlleles.get(0);
                }
            }

            // if still null, it doesn't exist...so create it
            if (canonicalAllele == null) {
                canonicalAllele = new CanonicalAllele();
                daoBean.getCanonicalAlleleDAO().save(canonicalAllele);
                canonicalAllele.getLocatedVariants().add(locatedVariant);
                if (liftOverLocatedVariant != null) {
                    canonicalAllele.getLocatedVariants().add(liftOverLocatedVariant);
                }
                daoBean.getCanonicalAlleleDAO().save(canonicalAllele);
            } else {

                if (!canonicalAllele.getLocatedVariants().contains(locatedVariant)) {
                    canonicalAllele.getLocatedVariants().add(locatedVariant);
                }

                if (liftOverLocatedVariant != null && !canonicalAllele.getLocatedVariants().contains(liftOverLocatedVariant)) {
                    canonicalAllele.getLocatedVariants().add(liftOverLocatedVariant);
                }
                daoBean.getCanonicalAlleleDAO().save(canonicalAllele);
            }
        } catch (CANVASDAOException e) {
            throw new BinningException(e);
        }

    }

    private void createAssmeblyLocatedVariantQC(String sampleName, VariantContext variantContext, LocatedVariant locatedVariant,
            Assembly assembly) throws BinningException {

        try {

            CommonInfo commonInfo = variantContext.getCommonInfo();

            Double qualityByDepth = commonInfo.hasAttribute("QD") ? Double.valueOf(commonInfo.getAttribute("QD").toString()) : null;

            Double readPosRankSum = commonInfo.hasAttribute("ReadPosRankSum")
                    ? Double.valueOf(commonInfo.getAttribute("ReadPosRankSum").toString()) : null;

            Integer homopolymerRun = commonInfo.hasAttribute("HRun") ? Integer.valueOf(commonInfo.getAttribute("HRun").toString()) : null;

            Double dels = commonInfo.hasAttribute("Dels") ? Double.valueOf(commonInfo.getAttribute("Dels").toString()) : null;

            Double fs = commonInfo.hasAttribute("FS") ? Double.valueOf(commonInfo.getAttribute("FS").toString()) : null;

            Genotype genotype = variantContext.getGenotype(sampleName);
            Integer dp = genotype.getDP();
            int[] ad = genotype.getAD();

            Integer refDepth = null;
            Integer altDepth = null;

            if (ad != null && ad.length == 2) {
                refDepth = ad[0];
                altDepth = ad[1];
            }

            AssemblyLocatedVariantPK alKey = new AssemblyLocatedVariantPK(assembly.getId(), locatedVariant.getId());
            AssemblyLocatedVariant alv = new AssemblyLocatedVariant(alKey);
            alv.setAssembly(assembly);
            alv.setLocatedVariant(locatedVariant);
            if (genotype.hasGQ()) {
                alv.setGenotypeQuality(Double.valueOf(genotype.getGQ() >= 0 ? genotype.getGQ() : -1));
            }
            alv.setHomozygous(genotype.isHom());
            daoBean.getAssemblyLocatedVariantDAO().save(alv);
            logger.info(alv.toString());

            // this is retarded...following database constraints
            if (dp != null || qualityByDepth != null || readPosRankSum != null || dels != null || homopolymerRun != null || fs != null
                    || refDepth != null || altDepth != null) {

                AssemblyLocatedVariantQCPK alvKey = new AssemblyLocatedVariantQCPK(assembly.getId(), locatedVariant.getId());
                AssemblyLocatedVariantQC alvQC = new AssemblyLocatedVariantQC(alvKey);
                alvQC.setAssembly(assembly);
                alvQC.setLocatedVariant(locatedVariant);

                if (dp == null || (dp != null && dp >= 0)) {
                    alvQC.setDepth(dp);
                }
                if (qualityByDepth == null || (qualityByDepth != null && qualityByDepth >= 0)) {
                    alvQC.setQualityByDepth(qualityByDepth);
                }
                alvQC.setReadPosRankSum(readPosRankSum);
                if (dels == null || (dels != null && dels >= 0 && dels <= 1)) {
                    alvQC.setFracReadsWithDels(dels);
                }
                if (homopolymerRun == null || (homopolymerRun != null && homopolymerRun >= 0)) {
                    alvQC.setHomopolymerRun(homopolymerRun);
                }
                if (fs == null || (fs != null && fs >= 0)) {
                    alvQC.setStrandScore(fs);
                }
                if (refDepth == null || (refDepth != null && refDepth >= 0)) {
                    alvQC.setRefDepth(refDepth);
                }
                if (altDepth == null || (altDepth != null && altDepth >= 0)) {
                    alvQC.setAltDepth(altDepth);
                }
                daoBean.getAssemblyLocatedVariantQCDAO().save(alvQC);
                logger.info(alvQC.toString());
            }
        } catch (NumberFormatException | CANVASDAOException e) {
            throw new BinningException(e);
        }

    }

    private Assembly buildAssembly(String sampleName) throws BinningException {
        logger.debug("ENTERING buildAssembly(String)");
        Assembly assembly = null;

        try {
            Lab lab = null;
            List<Lab> foundLabs = daoBean.getLabDAO().findByName(getLabName());
            if (CollectionUtils.isEmpty(foundLabs)) {
                lab = new Lab(getLabName());
                daoBean.getLabDAO().save(lab);
            } else {
                lab = foundLabs.get(0);
            }
            logger.info(lab.toString());

            Project project = daoBean.getProjectDAO().findById(getStudyName());
            if (project == null) {
                project = new Project(getStudyName());
                project.setLab(lab);
                daoBean.getProjectDAO().save(project);
                project = daoBean.getProjectDAO().findById(getStudyName());
            }
            logger.info(project.toString());

            Sample sample = null;
            List<Sample> foundSamples = daoBean.getSampleDAO().findByNameAndProjectName(binningJob.getParticipant(), project.getId());
            if (CollectionUtils.isEmpty(foundSamples)) {
                sample = new Sample(binningJob.getParticipant());
                sample.setProject(project);
                sample.setId(daoBean.getSampleDAO().save(sample));
            } else {
                sample = foundSamples.get(0);
            }
            logger.info(sample.toString());

            // library
            Library library = null;
            List<Library> foundLibraries = daoBean.getLibraryDAO().findByNameAndSampleId(getLibraryName(), sample.getId());
            if (CollectionUtils.isEmpty(foundLibraries)) {
                library = new Library(getLibraryName());
                library.setSample(sample);
                library.setId(daoBean.getLibraryDAO().save(library));
            } else {
                library = foundLibraries.get(0);
            }
            logger.info(library.toString());

            List<Assembly> foundAssemblies = daoBean.getAssemblyDAO().findByLibraryId(library.getId());
            if (CollectionUtils.isEmpty(foundAssemblies)) {

                VariantSet variantSet = new VariantSet();
                variantSet.setGenomeRef(getGenomeRef());
                variantSet.setId(daoBean.getVariantSetDAO().save(variantSet));

                assembly = new Assembly();
                assembly.setLibrary(library);
                assembly.setVariantSet(variantSet);
                assembly.setId(daoBean.getAssemblyDAO().save(assembly));

                binningJob.setAssembly(assembly);
                daoBean.getDiagnosticBinningJobDAO().save(binningJob);

            } else {
                assembly = binningJob.getAssembly();
                // delete asm loc var & qc instances
                logger.info("deleting AssemblyLocatedVariant instances");
                daoBean.getAssemblyLocatedVariantDAO().deleteByAssemblyId(binningJob.getAssembly().getId());
                logger.info("deleting AssemblyLocatedVariantQC instances");
                daoBean.getAssemblyLocatedVariantQCDAO().deleteByAssemblyId(binningJob.getAssembly().getId());
            }
        } catch (CANVASDAOException e) {
            e.printStackTrace();
        }
        return assembly;
    }

    public CANVASDAOBeanService getDaoBean() {
        return daoBean;
    }

    public void setDaoBean(CANVASDAOBeanService daoBean) {
        this.daoBean = daoBean;
    }

    public DiagnosticBinningJob getBinningJob() {
        return binningJob;
    }

    public void setBinningJob(DiagnosticBinningJob binningJob) {
        this.binningJob = binningJob;
    }

}
