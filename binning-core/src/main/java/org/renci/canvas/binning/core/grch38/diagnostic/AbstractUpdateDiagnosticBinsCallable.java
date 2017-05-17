package org.renci.canvas.binning.core.grch38.diagnostic;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.apache.commons.collections.CollectionUtils;
import org.renci.canvas.binning.core.BinningException;
import org.renci.canvas.binning.core.grch38.BinResultsFinalDiagnosticFactory;
import org.renci.canvas.binning.core.grch38.VariantsFactory;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.BinResultsFinalDiagnostic;
import org.renci.canvas.dao.clinbin.model.DiagnosticBinningJob;
import org.renci.canvas.dao.refseq.model.Variants_80_4;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public abstract class AbstractUpdateDiagnosticBinsCallable implements Callable<Void> {

    private static final Logger logger = LoggerFactory.getLogger(AbstractUpdateDiagnosticBinsCallable.class);

    private BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory;

    private CANVASDAOBeanService daoBean;

    private DiagnosticBinningJob binningJob;

    public AbstractUpdateDiagnosticBinsCallable(CANVASDAOBeanService daoBean, DiagnosticBinningJob binningJob) {
        super();
        this.daoBean = daoBean;
        this.binningJob = binningJob;
        binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
    }

    @Override
    public Void call() throws BinningException {
        logger.debug("ENTERING call()");

        try {

            logger.info("Deleting BinResultsFinalDiagnostic instances by assembly id = {}", binningJob.getAssembly().getId());
            daoBean.getBinResultsFinalDiagnosticDAO().deleteByAssemblyId(binningJob.getAssembly().getId());

            List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO().findByAssemblyId(binningJob.getAssembly().getId());

            List<Variants_80_4> variants = new ArrayList<>();

            if (CollectionUtils.isNotEmpty(locatedVariantList)) {
                logger.info(String.format("locatedVariantList.size(): %d", locatedVariantList.size()));
                for (LocatedVariant locatedVariant : locatedVariantList) {
                    List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                    if (CollectionUtils.isNotEmpty(foundVariants)) {
                        variants.addAll(foundVariants);
                    }
                }
            }

            if (CollectionUtils.isNotEmpty(variants)) {

                logger.info(String.format("variants.size(): %d", variants.size()));

                try {

                    // hgmd known pathogenic...disease class 1
                    List<BinResultsFinalDiagnostic> knownPathogenic = binResultsFinalDiagnosticFactory.findHGMDKnownPathogenic(binningJob,
                            variants);

                    if (CollectionUtils.isNotEmpty(knownPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : knownPathogenic) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 1);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar known pathogenic...disease class 1
                    knownPathogenic = binResultsFinalDiagnosticFactory.findClinVarKnownPathogenic(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(knownPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : knownPathogenic) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 1);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }

                try {
                    // hgmd likely pathogenic...disease class 2
                    List<BinResultsFinalDiagnostic> likelyPathogenic = binResultsFinalDiagnosticFactory.findHGMDLikelyPathogenic(binningJob,
                            variants);
                    if (CollectionUtils.isNotEmpty(likelyPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : likelyPathogenic) {
                            logger.info(binResultsFinalDiagnostic.getId().toString());
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 2);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar likely pathogenic...disease class 2
                    likelyPathogenic = binResultsFinalDiagnosticFactory.findHGMDLikelyPathogenic(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(likelyPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : likelyPathogenic) {
                            logger.info(binResultsFinalDiagnostic.getId().toString());
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 2);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }

                try {
                    // hgmd possibly pathogenic...disease class 3
                    List<BinResultsFinalDiagnostic> possiblyPathogenic = binResultsFinalDiagnosticFactory
                            .findHGMDPossiblyPathogenic(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(possiblyPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : possiblyPathogenic) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 3);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar possibly pathogenic...disease class 3
                    possiblyPathogenic = binResultsFinalDiagnosticFactory.findClinVarPossiblyPathogenic(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(possiblyPathogenic)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : possiblyPathogenic) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 3);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }

                try {
                    // hgmd uncertain significance...disease class 4
                    List<BinResultsFinalDiagnostic> uncertainSignificance = binResultsFinalDiagnosticFactory
                            .findHGMDUncertainSignificance(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(uncertainSignificance)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : uncertainSignificance) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 4);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar uncertain significance...disease class 4
                    uncertainSignificance = binResultsFinalDiagnosticFactory.findClinVarUncertainSignificance(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(uncertainSignificance)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : uncertainSignificance) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 4);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }

                try {
                    // hgmd likely benign...disease class 5
                    List<BinResultsFinalDiagnostic> likelyBenign = binResultsFinalDiagnosticFactory.findHGMDLikelyBenign(binningJob,
                            variants);
                    if (CollectionUtils.isNotEmpty(likelyBenign)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : likelyBenign) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 5);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar likely benign...disease class 5
                    likelyBenign = binResultsFinalDiagnosticFactory.findHGMDLikelyBenign(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(likelyBenign)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : likelyBenign) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 5);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }

                try {
                    // hgmd almost certainly benign...disease class 6
                    List<BinResultsFinalDiagnostic> almostCertainlyBenign = binResultsFinalDiagnosticFactory
                            .findHGMDAlmostCertainlyBenign(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(almostCertainlyBenign)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : almostCertainlyBenign) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndHGMDDiseaseClassId(binResultsFinalDiagnostic.getId(), 6);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }

                    // clinvar almost certainly benign...disease class 6
                    almostCertainlyBenign = binResultsFinalDiagnosticFactory.findClinVarAlmostCertainlyBenign(binningJob, variants);
                    if (CollectionUtils.isNotEmpty(almostCertainlyBenign)) {
                        for (BinResultsFinalDiagnostic binResultsFinalDiagnostic : almostCertainlyBenign) {
                            List<BinResultsFinalDiagnostic> foundBinResultsFinalDiagnostics = daoBean.getBinResultsFinalDiagnosticDAO()
                                    .findByKeyAndClinVarDiseaseClassId(binResultsFinalDiagnostic.getId(), 6);
                            if (CollectionUtils.isEmpty(foundBinResultsFinalDiagnostics)) {
                                logger.info("saving BinResultsFinalDiagnostic: {}", binResultsFinalDiagnostic.toString());
                                // daoBean.getBinResultsFinalDiagnosticDAO().save(binResultsFinalDiagnostic);
                            }
                        }
                    }
                } catch (CANVASDAOException e) {
                    logger.error(e.getMessage(), e);
                }
                // });
                //
                // es.shutdown();
                // es.awaitTermination(2L, TimeUnit.HOURS);
            }

        } catch (Exception e) {
            logger.error(e.getMessage(), e);
            throw new BinningException(e);
        }

        return null;
    }

}
