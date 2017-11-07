package org.renci.canvas.binning.core.grch38.diagnostic;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.apache.commons.collections.CollectionUtils;
import org.renci.canvas.binning.core.BinningException;
import org.renci.canvas.binning.core.grch38.VariantsFactory;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.DiagnosticBinningJob;
import org.renci.canvas.dao.clinbin.model.DiagnosticResultVersion;
import org.renci.canvas.dao.ref.model.GenomeRef;
import org.renci.canvas.dao.refseq.Variants_80_4_DAO;
import org.renci.canvas.dao.refseq.model.TranscriptMaps;
import org.renci.canvas.dao.refseq.model.TranscriptMapsExons;
import org.renci.canvas.dao.refseq.model.Variants_80_4;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class AbstractAnnotateVariantsCallable implements Callable<Void> {

    private static final Logger logger = LoggerFactory.getLogger(AbstractAnnotateVariantsCallable.class);

    private CANVASDAOBeanService daoBean;

    private DiagnosticBinningJob binningJob;

    public AbstractAnnotateVariantsCallable(CANVASDAOBeanService daoBean, DiagnosticBinningJob binningJob) {
        super();
        this.daoBean = daoBean;
        this.binningJob = binningJob;
    }

    @Override
    public Void call() throws BinningException {
        logger.debug("ENTERING run()");

        try {

            DiagnosticResultVersion diagnosticResultVersion = binningJob.getDiagnosticResultVersion();
            logger.info(diagnosticResultVersion.toString());

            String refseqVersion = diagnosticResultVersion.getRefseqVersion().toString();

            GenomeRef genomeRef = diagnosticResultVersion.getGenomeRef();
            logger.info(genomeRef.toString());

            List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO().findByAssemblyId(binningJob.getAssembly().getId());

            if (CollectionUtils.isNotEmpty(locatedVariantList)) {
                logger.info(String.format("locatedVariantList.size(): %d", locatedVariantList.size()));

                ExecutorService es = Executors.newFixedThreadPool(4);
                for (LocatedVariant locatedVariant : locatedVariantList) {

                    es.submit(() -> updateVariant(locatedVariant, refseqVersion, genomeRef.getId(), daoBean));

                }

                es.shutdown();
                if (!es.awaitTermination(1L, TimeUnit.DAYS)) {
                    es.shutdownNow();
                }

            }

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            throw new BinningException(e);
        }

        return null;
    }

    private void updateVariant(LocatedVariant locatedVariant, String refseqVersion, Integer genomeRefId, CANVASDAOBeanService daoBean) {
        Variants_80_4_DAO dao = daoBean.getVariants_80_4_DAO();
        try {
            dao.deleteByLocatedVariantId(locatedVariant.getId());
        } catch (CANVASDAOException e) {
            logger.error(e.getMessage(), e);
        }
        Collection<Variants_80_4> variants = annotateVariant(locatedVariant, refseqVersion, genomeRefId, daoBean);
        for (Variants_80_4 v : variants) {
            try {
                dao.save(v);
            } catch (CANVASDAOException e) {
                logger.error(e.getMessage(), e);
            }
        }

    }

    public static Set<Variants_80_4> annotateVariant(LocatedVariant locatedVariant, String refseqVersion, Integer genomeRefId, CANVASDAOBeanService daoBean) {

        logger.info(locatedVariant.toString());
        VariantsFactory variantsFactory = VariantsFactory.getInstance(daoBean);
        Set<Variants_80_4> variants = new HashSet<>();

        try {

            final List<TranscriptMaps> transcriptMapsList = daoBean.getTranscriptMapsDAO()
                    .findByGenomeRefIdAndRefSeqVersionAndGenomeRefSeqAccessionAndInExonRange(genomeRefId,
                            refseqVersion, locatedVariant.getGenomeRefSeq().getId(), locatedVariant.getPosition());

            // FIXME - this still only finds transcripts where the variant at least starts or ends-- not variants that span the whole transcript...
            // FIXME - would be better to have a proper overlap query
            transcriptMapsList.addAll(daoBean.getTranscriptMapsDAO()
                    .findByGenomeRefIdAndRefSeqVersionAndGenomeRefSeqAccessionAndInExonRange(genomeRefId,
                            refseqVersion, locatedVariant.getGenomeRefSeq().getId(), locatedVariant.getEndPosition()-1));

            // now we need to use distinct() to eliminate duplicates (e.g., variant starts and ends within transcript)
            List<TranscriptMaps> distinctTranscriptMapsList = transcriptMapsList.stream()
                    .map(a -> a.getTranscript().getId())
                    .distinct().map(a -> transcriptMapsList.parallelStream()
                            .filter(b -> b.getTranscript().getId().equals(a)).findAny().get())
                    .collect(Collectors.toList());

            if (CollectionUtils.isNotEmpty(distinctTranscriptMapsList)) {

                distinctTranscriptMapsList.sort((a, b) -> b.getTranscript().getId().compareTo(a.getTranscript().getId()));

                for (TranscriptMaps tMap : distinctTranscriptMapsList) {

                    logger.info(tMap.toString());

                    // I think this is only used to provide mapnum and nummaps
                    List<TranscriptMaps> mapsList = daoBean.getTranscriptMapsDAO()
                            .findByGenomeRefIdAndRefSeqVersionAndTranscriptId(genomeRefId, refseqVersion,
                                    tMap.getTranscript().getId());
                    mapsList.sort((a, b) -> b.getId().compareTo(a.getId()));

                    List<TranscriptMapsExons> transcriptMapsExonsList = daoBean.getTranscriptMapsExonsDAO()
                            .findByTranscriptMapsId(tMap.getId());

                    long varEnd = locatedVariant.getPosition()+locatedVariant.getRef().length()-1; // FIXME - can we just use getEndPosition()-1?
                    boolean foundExonOverlap = false;
                    for (TranscriptMapsExons exon : transcriptMapsExonsList) {
                        // FIXME - might be able to `break` depending on strandedness..
                        long exonMin = Math.min(exon.getContigStart(), exon.getContigEnd());
                        long exonMax = Math.max(exon.getContigStart(), exon.getContigEnd());
                        if (varEnd < exonMin) continue;
                        if (locatedVariant.getPosition() > exonMax) continue;
                        // we have overlap...
                        foundExonOverlap = true;
                        if (locatedVariant.getPosition() >= exonMin
                                && varEnd <= exonMax) {
                            // completely contained
                            // FIXME - have to shallow copy transcriptsMapsExonsList because something in there reorders it...
                            variants.add(variantsFactory.createExonicVariant(locatedVariant, mapsList,
                                    new ArrayList<TranscriptMapsExons>(transcriptMapsExonsList), exon));
                        } else {
                            variants.add(variantsFactory.createBorderCrossingVariant(locatedVariant, tMap, mapsList, transcriptMapsExonsList, exon));
                        }
                    }

                    if (!foundExonOverlap) {
                        variants.add(variantsFactory.createIntronicVariant(locatedVariant, mapsList, tMap, transcriptMapsExonsList));
                    }
                }
            }

            if (CollectionUtils.isEmpty(variants)) {
                // not found in or across any transcript, must be intergenic
                Variants_80_4 variant = variantsFactory.createIntergenicVariant(locatedVariant);
                variants.add(variant);
            }

            for (Variants_80_4 variant : variants) {
                logger.info(variant.toString());
            }

        } catch (CANVASDAOException | BinningException e) {
            logger.error(e.getMessage(), e);
        }
        return variants;

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
