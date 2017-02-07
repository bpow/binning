package org.renci.binning.core.incidental;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

import org.apache.commons.collections.CollectionUtils;
import org.apache.commons.lang3.Range;
import org.renci.binning.core.BinningException;
import org.renci.binning.core.grch37.VariantsFactory;
import org.renci.binning.dao.BinningDAOBeanService;
import org.renci.binning.dao.clinbin.model.IncidentalBinningJob;
import org.renci.binning.dao.clinbin.model.IncidentalResultVersionX;
import org.renci.binning.dao.ref.model.GenomeRef;
import org.renci.binning.dao.refseq.model.TranscriptMaps;
import org.renci.binning.dao.refseq.model.TranscriptMapsExons;
import org.renci.binning.dao.refseq.model.Variants_61_2;
import org.renci.binning.dao.var.model.LocatedVariant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class AbstractAnnotateVariantsCallable implements Callable<List<Variants_61_2>> {

    private static final Logger logger = LoggerFactory.getLogger(AbstractAnnotateVariantsCallable.class);

    private BinningDAOBeanService daoBean;

    private IncidentalBinningJob binningJob;

    public AbstractAnnotateVariantsCallable(BinningDAOBeanService daoBean, IncidentalBinningJob binningJob) {
        super();
        this.daoBean = daoBean;
        this.binningJob = binningJob;
    }

    @Override
    public List<Variants_61_2> call() throws BinningException {
        logger.debug("ENTERING run()");

        List<Variants_61_2> variants = new ArrayList<>();

        try {

            IncidentalResultVersionX diagnosticResultVersion = daoBean.getIncidentalResultVersionXDAO()
                    .findById(binningJob.getListVersion());
            logger.info(diagnosticResultVersion.toString());

            String refseqVersion = diagnosticResultVersion.getRefseqVersion().toString();

            GenomeRef genomeRef = diagnosticResultVersion.getGenomeRef();
            logger.info(genomeRef.toString());

            List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO().findByAssemblyId(binningJob.getAssembly().getId());

            if (CollectionUtils.isNotEmpty(locatedVariantList)) {
                logger.info(String.format("locatedVariantList.size(): %d", locatedVariantList.size()));
                locatedVariantList.sort((a, b) -> {
                    int ret = a.getGenomeRefSeq().getVerAccession().compareTo(b.getGenomeRefSeq().getVerAccession());
                    if (ret == 0) {
                        ret = a.getPosition().compareTo(b.getPosition());
                    }
                    return ret;
                });
                for (LocatedVariant locatedVariant : locatedVariantList) {
                    logger.info(locatedVariant.toString());

                    List<TranscriptMaps> transcriptMapsList = daoBean.getTranscriptMapsDAO()
                            .findByGenomeRefIdAndRefSeqVersionAndGenomeRefSeqAccessionAndInExonRange(genomeRef.getId(), refseqVersion,
                                    locatedVariant.getGenomeRefSeq().getVerAccession(), locatedVariant.getPosition());

                    if (CollectionUtils.isNotEmpty(transcriptMapsList)) {

                        Map<String, List<TranscriptMaps>> transcriptMap = new HashMap<String, List<TranscriptMaps>>();
                        for (TranscriptMaps tMap : transcriptMapsList) {
                            if (!transcriptMap.containsKey(tMap.getTranscript().getVersionId())) {
                                transcriptMap.put(tMap.getTranscript().getVersionId(), new ArrayList<TranscriptMaps>());
                            }
                            transcriptMap.get(tMap.getTranscript().getVersionId()).add(tMap);
                        }
                        TranscriptMaps remove = null;
                        for (String key : transcriptMap.keySet()) {
                            if (transcriptMap.get(key).size() > 1) {
                                transcriptMap.get(key).sort((a, b) -> a.getMapCount().compareTo(b.getMapCount()));
                                remove = transcriptMap.get(key).get(0);
                                logger.info("removing duplicate: {}", remove.toString());
                                transcriptMapsList.remove(remove);
                            }
                        }
                        transcriptMapsList.sort((a, b) -> b.getTranscript().getVersionId().compareTo(a.getTranscript().getVersionId()));

                        // handling non boundary crossing variants (intron/exon/utr*)
                        logger.info("transcriptMapsList.size(): {}", transcriptMapsList.size());

                        for (TranscriptMaps tMap : transcriptMapsList) {

                            logger.info(tMap.toString());
                            List<TranscriptMapsExons> transcriptMapsExonsList = daoBean.getTranscriptMapsExonsDAO()
                                    .findByTranscriptMapsId(tMap.getId());

                            TranscriptMapsExons transcriptMapsExons = null;
                            if (CollectionUtils.isNotEmpty(transcriptMapsExonsList)) {
                                for (TranscriptMapsExons exon : transcriptMapsExonsList) {
                                    Range<Integer> contigRange = Range.between(exon.getContigStart(), exon.getContigEnd());
                                    if (contigRange.contains(locatedVariant.getPosition())) {
                                        transcriptMapsExons = exon;
                                        logger.info(transcriptMapsExons.toString());
                                        break;
                                    }
                                }
                            }

                            List<TranscriptMaps> mapsList = daoBean.getTranscriptMapsDAO().findByGenomeRefIdAndRefSeqVersionAndTranscriptId(
                                    genomeRef.getId(), refseqVersion, tMap.getTranscript().getVersionId());

                            Variants_61_2 variant = null;

                            if (transcriptMapsExons == null) {
                                variant = VariantsFactory.createIntronicVariant(daoBean, refseqVersion, locatedVariant, mapsList, tMap,
                                        transcriptMapsExonsList);
                            } else {

                                if (!"snp".equals(locatedVariant.getVariantType().getName())
                                        && ((transcriptMapsExons.getContigEnd().equals(locatedVariant.getPosition())
                                                && "-".equals(tMap.getStrand()))
                                                || (transcriptMapsExons.getContigStart().equals(locatedVariant.getPosition())
                                                        && "+".equals(tMap.getStrand())))) {
                                    variant = VariantsFactory.createBorderCrossingVariant(daoBean, refseqVersion, locatedVariant, tMap,
                                            mapsList, transcriptMapsExonsList, transcriptMapsExons);
                                } else {
                                    variant = VariantsFactory.createExonicVariant(daoBean, refseqVersion, locatedVariant, mapsList,
                                            transcriptMapsExonsList, transcriptMapsExons);
                                }

                            }
                            variants.add(variant);
                        }

                    } else {

                        // try searching by adjusting for length of locatedVariant.getSeq()...could be intron/exon
                        // boundary crossing
                        transcriptMapsList = daoBean.getTranscriptMapsDAO()
                                .findByGenomeRefIdAndRefSeqVersionAndGenomeRefSeqAccessionAndInExonRange(genomeRef.getId(), refseqVersion,
                                        locatedVariant.getGenomeRefSeq().getVerAccession(),
                                        locatedVariant.getPosition() + locatedVariant.getRef().length() - 1);

                        if (CollectionUtils.isNotEmpty(transcriptMapsList)) {

                            Map<String, List<TranscriptMaps>> transcriptMap = new HashMap<String, List<TranscriptMaps>>();
                            for (TranscriptMaps tMap : transcriptMapsList) {
                                if (!transcriptMap.containsKey(tMap.getTranscript().getVersionId())) {
                                    transcriptMap.put(tMap.getTranscript().getVersionId(), new ArrayList<TranscriptMaps>());
                                }
                                transcriptMap.get(tMap.getTranscript().getVersionId()).add(tMap);
                            }
                            TranscriptMaps remove = null;
                            for (String key : transcriptMap.keySet()) {
                                if (transcriptMap.get(key).size() > 1) {
                                    transcriptMap.get(key).sort((a, b) -> a.getMapCount().compareTo(b.getMapCount()));
                                    remove = transcriptMap.get(key).get(0);
                                    logger.info("removing duplicate: {}", remove.toString());
                                    transcriptMapsList.remove(remove);
                                }
                            }
                            transcriptMapsList.sort((a, b) -> b.getTranscript().getVersionId().compareTo(a.getTranscript().getVersionId()));

                            for (TranscriptMaps tMap : transcriptMapsList) {
                                logger.info(tMap.toString());

                                List<TranscriptMapsExons> transcriptMapsExonsList = daoBean.getTranscriptMapsExonsDAO()
                                        .findByTranscriptMapsId(tMap.getId());

                                TranscriptMapsExons transcriptMapsExons = null;
                                if (CollectionUtils.isNotEmpty(transcriptMapsExonsList)) {
                                    for (TranscriptMapsExons exon : transcriptMapsExonsList) {
                                        Range<Integer> contigRange = Range.between(exon.getContigStart(), exon.getContigEnd());
                                        if (contigRange.contains(locatedVariant.getPosition() + locatedVariant.getRef().length() - 1)) {
                                            transcriptMapsExons = exon;
                                            logger.info(transcriptMapsExons.toString());
                                            break;
                                        }
                                    }
                                }

                                List<TranscriptMaps> mapsList = daoBean.getTranscriptMapsDAO()
                                        .findByGenomeRefIdAndRefSeqVersionAndTranscriptId(genomeRef.getId(), refseqVersion,
                                                tMap.getTranscript().getVersionId());

                                Variants_61_2 variant = VariantsFactory.createBorderCrossingVariant(daoBean, refseqVersion,
                                        locatedVariant, tMap, mapsList, transcriptMapsExonsList, transcriptMapsExons);
                                variants.add(variant);
                            }
                        } else {
                            transcriptMapsList = daoBean.getTranscriptMapsDAO()
                                    .findByGenomeRefIdAndRefSeqVersionAndGenomeRefSeqAccessionAndInExonRange(genomeRef.getId(),
                                            refseqVersion, locatedVariant.getGenomeRefSeq().getVerAccession(),
                                            locatedVariant.getPosition() - locatedVariant.getRef().length());

                            if (CollectionUtils.isNotEmpty(transcriptMapsList)) {
                                Map<String, List<TranscriptMaps>> transcriptMap = new HashMap<String, List<TranscriptMaps>>();
                                for (TranscriptMaps tMap : transcriptMapsList) {
                                    if (!transcriptMap.containsKey(tMap.getTranscript().getVersionId())) {
                                        transcriptMap.put(tMap.getTranscript().getVersionId(), new ArrayList<TranscriptMaps>());
                                    }
                                    transcriptMap.get(tMap.getTranscript().getVersionId()).add(tMap);
                                }
                                TranscriptMaps remove = null;
                                for (String key : transcriptMap.keySet()) {
                                    if (transcriptMap.get(key).size() > 1) {
                                        transcriptMap.get(key).sort((a, b) -> a.getMapCount().compareTo(b.getMapCount()));
                                        remove = transcriptMap.get(key).get(0);
                                        logger.info("removing duplicate: {}", remove.toString());
                                        transcriptMapsList.remove(remove);
                                    }
                                }
                                transcriptMapsList
                                        .sort((a, b) -> b.getTranscript().getVersionId().compareTo(a.getTranscript().getVersionId()));

                                for (TranscriptMaps tMap : transcriptMapsList) {
                                    logger.info(tMap.toString());

                                    List<TranscriptMapsExons> transcriptMapsExonsList = daoBean.getTranscriptMapsExonsDAO()
                                            .findByTranscriptMapsId(tMap.getId());

                                    TranscriptMapsExons transcriptMapsExons = null;
                                    if (CollectionUtils.isNotEmpty(transcriptMapsExonsList)) {
                                        for (TranscriptMapsExons exon : transcriptMapsExonsList) {
                                            Range<Integer> contigRange = Range.between(exon.getContigStart(), exon.getContigEnd());
                                            if (contigRange.contains(locatedVariant.getPosition() - locatedVariant.getRef().length())) {
                                                transcriptMapsExons = exon;
                                                logger.info(transcriptMapsExons.toString());
                                                break;
                                            }
                                        }
                                    }

                                    List<TranscriptMaps> mapsList = daoBean.getTranscriptMapsDAO()
                                            .findByGenomeRefIdAndRefSeqVersionAndTranscriptId(genomeRef.getId(), refseqVersion,
                                                    tMap.getTranscript().getVersionId());

                                    Variants_61_2 variant = VariantsFactory.createBorderCrossingVariant(daoBean, refseqVersion,
                                            locatedVariant, tMap, mapsList, transcriptMapsExonsList, transcriptMapsExons);
                                    variants.add(variant);
                                }
                            }
                        }

                    }

                }

            }

            if (CollectionUtils.isNotEmpty(variants)) {

                List<Variants_61_2> toRemove = new ArrayList<>();

                for (Variants_61_2 variant : variants) {
                    Variants_61_2 foundVariant = daoBean.getVariants_61_2_DAO().findById(variant.getKey());
                    if (foundVariant != null) {
                        toRemove.add(variant);
                    }
                }

                for (Variants_61_2 variant : toRemove) {
                    variants.remove(variant);
                }

            }

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            throw new BinningException(e);
        }

        return variants;
    }

    public BinningDAOBeanService getDaoBean() {
        return daoBean;
    }

    public void setDaoBean(BinningDAOBeanService daoBean) {
        this.daoBean = daoBean;
    }

    public IncidentalBinningJob getBinningJob() {
        return binningJob;
    }

    public void setBinningJob(IncidentalBinningJob binningJob) {
        this.binningJob = binningJob;
    }

}