package org.renci.binning.core;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import javax.persistence.EntityManager;
import javax.persistence.EntityManagerFactory;
import javax.persistence.Persistence;

import org.apache.commons.collections.CollectionUtils;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.renci.canvas.binning.core.grch38.BinResultsFinalDiagnosticFactory;
import org.renci.canvas.dao.CANVASDAOException;
import org.renci.canvas.dao.clinbin.model.BinResultsFinalDiagnostic;
import org.renci.canvas.dao.clinbin.model.DiagnosticBinningJob;
import org.renci.canvas.dao.dbsnp.model.SNPMappingAgg;
import org.renci.canvas.dao.jpa.CANVASDAOBeanServiceImpl;
import org.renci.canvas.dao.jpa.annotation.AnnotationGeneExternalIdDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.DiagnosticBinningJobDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.DiagnosticGeneDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.DiagnosticResultVersionDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.DiseaseClassDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.MaxFrequencyDAOImpl;
import org.renci.canvas.dao.jpa.clinbin.NCGenesFrequenciesDAOImpl;
import org.renci.canvas.dao.jpa.clinvar.ReferenceClinicalAssertionDAOImpl;
import org.renci.canvas.dao.jpa.dbsnp.SNPMappingAggDAOImpl;
import org.renci.canvas.dao.jpa.hgnc.HGNCGeneDAOImpl;
import org.renci.canvas.dao.jpa.refseq.FeatureDAOImpl;
import org.renci.canvas.dao.jpa.refseq.LocationTypeDAOImpl;
import org.renci.canvas.dao.jpa.refseq.RefSeqCodingSequenceDAOImpl;
import org.renci.canvas.dao.jpa.refseq.RefSeqGeneDAOImpl;
import org.renci.canvas.dao.jpa.refseq.RegionGroupRegionDAOImpl;
import org.renci.canvas.dao.jpa.refseq.TranscriptMapsDAOImpl;
import org.renci.canvas.dao.jpa.refseq.TranscriptMapsExonsDAOImpl;
import org.renci.canvas.dao.jpa.refseq.VariantEffectDAOImpl;
import org.renci.canvas.dao.jpa.refseq.Variants_80_4_DAOImpl;
import org.renci.canvas.dao.jpa.var.CanonicalAlleleDAOImpl;
import org.renci.canvas.dao.jpa.var.LocatedVariantDAOImpl;
import org.renci.canvas.dao.refseq.model.Variants_80_4;
import org.renci.canvas.dao.var.model.LocatedVariant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class UpdateDiagnosticBinsCallable38Test {

    private static final Logger logger = LoggerFactory.getLogger(UpdateDiagnosticBinsCallable38Test.class);

    private static EntityManagerFactory emf;

    private static EntityManager em;

    private static CANVASDAOBeanServiceImpl daoBean = new CANVASDAOBeanServiceImpl();

    public UpdateDiagnosticBinsCallable38Test() {
        super();
    }

    @BeforeClass
    public static void setup() {
        emf = Persistence.createEntityManagerFactory("canvas_test", null);
        em = emf.createEntityManager();

        TranscriptMapsDAOImpl transcriptMapsDAO = new TranscriptMapsDAOImpl();
        transcriptMapsDAO.setEntityManager(em);
        daoBean.setTranscriptMapsDAO(transcriptMapsDAO);

        TranscriptMapsExonsDAOImpl transcriptMapsExonsDAO = new TranscriptMapsExonsDAOImpl();
        transcriptMapsExonsDAO.setEntityManager(em);
        daoBean.setTranscriptMapsExonsDAO(transcriptMapsExonsDAO);

        LocationTypeDAOImpl locationTypeDAO = new LocationTypeDAOImpl();
        locationTypeDAO.setEntityManager(em);
        daoBean.setLocationTypeDAO(locationTypeDAO);

        VariantEffectDAOImpl variantEffectDAO = new VariantEffectDAOImpl();
        variantEffectDAO.setEntityManager(em);
        daoBean.setVariantEffectDAO(variantEffectDAO);

        RefSeqGeneDAOImpl refSeqGeneDAO = new RefSeqGeneDAOImpl();
        refSeqGeneDAO.setEntityManager(em);
        daoBean.setRefSeqGeneDAO(refSeqGeneDAO);

        HGNCGeneDAOImpl hgncGeneDAO = new HGNCGeneDAOImpl();
        hgncGeneDAO.setEntityManager(em);
        daoBean.setHGNCGeneDAO(hgncGeneDAO);

        AnnotationGeneExternalIdDAOImpl annotationGeneExternalIdDAO = new AnnotationGeneExternalIdDAOImpl();
        annotationGeneExternalIdDAO.setEntityManager(em);
        daoBean.setAnnotationGeneExternalIdDAO(annotationGeneExternalIdDAO);

        RegionGroupRegionDAOImpl regionGroupRegionDAO = new RegionGroupRegionDAOImpl();
        regionGroupRegionDAO.setEntityManager(em);
        daoBean.setRegionGroupRegionDAO(regionGroupRegionDAO);

        RefSeqCodingSequenceDAOImpl refSeqCodingSequenceDAO = new RefSeqCodingSequenceDAOImpl();
        refSeqCodingSequenceDAO.setEntityManager(em);
        daoBean.setRefSeqCodingSequenceDAO(refSeqCodingSequenceDAO);

        FeatureDAOImpl featureDAO = new FeatureDAOImpl();
        featureDAO.setEntityManager(em);
        daoBean.setFeatureDAO(featureDAO);

        LocatedVariantDAOImpl locatedVariantDAO = new LocatedVariantDAOImpl();
        locatedVariantDAO.setEntityManager(em);
        daoBean.setLocatedVariantDAO(locatedVariantDAO);

        DiagnosticBinningJobDAOImpl diagnosticBinningJobDAO = new DiagnosticBinningJobDAOImpl();
        diagnosticBinningJobDAO.setEntityManager(em);
        daoBean.setDiagnosticBinningJobDAO(diagnosticBinningJobDAO);

        Variants_80_4_DAOImpl variants_80_4_DAO = new Variants_80_4_DAOImpl();
        variants_80_4_DAO.setEntityManager(em);
        daoBean.setVariants_80_4_DAO(variants_80_4_DAO);

        DiseaseClassDAOImpl diseaseClassDAO = new DiseaseClassDAOImpl();
        diseaseClassDAO.setEntityManager(em);
        daoBean.setDiseaseClassDAO(diseaseClassDAO);

        NCGenesFrequenciesDAOImpl ncgenesFrequenciesDAO = new NCGenesFrequenciesDAOImpl();
        ncgenesFrequenciesDAO.setEntityManager(em);
        daoBean.setNCGenesFrequenciesDAO(ncgenesFrequenciesDAO);

        DiagnosticResultVersionDAOImpl diagnosticResultVersionDAO = new DiagnosticResultVersionDAOImpl();
        diagnosticResultVersionDAO.setEntityManager(em);
        daoBean.setDiagnosticResultVersionDAO(diagnosticResultVersionDAO);

        ReferenceClinicalAssertionDAOImpl referenceClinicalAssertionDAO = new ReferenceClinicalAssertionDAOImpl();
        referenceClinicalAssertionDAO.setEntityManager(em);
        daoBean.setReferenceClinicalAssertionDAO(referenceClinicalAssertionDAO);

        CanonicalAlleleDAOImpl canonicalAlleleDAO = new CanonicalAlleleDAOImpl();
        canonicalAlleleDAO.setEntityManager(em);
        daoBean.setCanonicalAlleleDAO(canonicalAlleleDAO);

        MaxFrequencyDAOImpl maxFrequencyDAO = new MaxFrequencyDAOImpl();
        maxFrequencyDAO.setEntityManager(em);
        daoBean.setMaxFrequencyDAO(maxFrequencyDAO);

        DiagnosticGeneDAOImpl diagnosticGeneDAO = new DiagnosticGeneDAOImpl();
        diagnosticGeneDAO.setEntityManager(em);
        daoBean.setDiagnosticGeneDAO(diagnosticGeneDAO);

        SNPMappingAggDAOImpl snpMappingAggDAO = new SNPMappingAggDAOImpl();
        snpMappingAggDAO.setEntityManager(em);
        daoBean.setSNPMappingAggDAO(snpMappingAggDAO);

    }

    @AfterClass
    public static void tearDown() {
        em.close();
        emf.close();
    }

    @Test
    public void testKnownPathogenic() throws CANVASDAOException {

        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);

        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(5001);

        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                // assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                if (CollectionUtils.isNotEmpty(foundVariants)) {
                    results.addAll(binResultsFinalDiagnosticFactory.findClinVarKnownPathogenic(diagnosticBinningJob, foundVariants));
                }
            }
        }

        assertTrue(CollectionUtils.isNotEmpty(results));
        results.forEach(a -> logger.info(a.toString()));
        logger.info("results.size(): {}", results.size());
    }

    @Test
    public void testLikelyPathenogenic() throws CANVASDAOException {
        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(4218);
        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            List<Variants_80_4> variants = new ArrayList<>();
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                variants.addAll(foundVariants);
            }
            results.addAll(binResultsFinalDiagnosticFactory.findHGMDLikelyPathogenic(diagnosticBinningJob, variants));
        }
        assertTrue(CollectionUtils.isNotEmpty(results));
        results.forEach(a -> logger.info(a.toString()));
        logger.info("results.size(): {}", results.size());
        assertTrue(results.size() == 11);
    }

    @Test
    public void testPossiblyPathenogenic() throws CANVASDAOException {
        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(4207);
        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            List<Variants_80_4> variants = new ArrayList<>();
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                variants.addAll(foundVariants);
            }
            results.addAll(binResultsFinalDiagnosticFactory.findHGMDPossiblyPathogenic(diagnosticBinningJob, variants));
        }
        assertTrue(CollectionUtils.isNotEmpty(results));
        logger.info("results.size(): {}", results.size());
        assertTrue(results.size() == 54);
        results.forEach(a -> logger.info(a.toString()));
    }

    @Test
    public void testUncertainSignificance() throws CANVASDAOException {
        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(4207);
        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            List<Variants_80_4> variants = new ArrayList<>();
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                variants.addAll(foundVariants);
            }
            results.addAll(binResultsFinalDiagnosticFactory.findHGMDUncertainSignificance(diagnosticBinningJob, variants));
        }
        assertTrue(CollectionUtils.isNotEmpty(results));
        logger.info("results.size(): {}", results.size());
        assertTrue(results.size() == 151);
        results.forEach(a -> logger.info(a.toString()));
    }

    @Test
    public void testLikelyBenign() throws CANVASDAOException {
        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(4207);
        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            List<Variants_80_4> variants = new ArrayList<>();
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                variants.addAll(foundVariants);
            }
            results.addAll(binResultsFinalDiagnosticFactory.findHGMDLikelyBenign(diagnosticBinningJob, variants));
        }
        assertTrue(CollectionUtils.isNotEmpty(results));
        logger.info("results.size(): {}", results.size());
        assertTrue(results.size() == 1);
        results.forEach(a -> logger.info(a.toString()));
    }

    @Test
    public void testAlmostCertainlyBenign() throws CANVASDAOException {
        BinResultsFinalDiagnosticFactory binResultsFinalDiagnosticFactory = BinResultsFinalDiagnosticFactory.getInstance(daoBean);
        DiagnosticBinningJob diagnosticBinningJob = daoBean.getDiagnosticBinningJobDAO().findById(4207);
        List<BinResultsFinalDiagnostic> results = new ArrayList<>();
        List<LocatedVariant> locatedVariantList = daoBean.getLocatedVariantDAO()
                .findByAssemblyId(diagnosticBinningJob.getAssembly().getId());

        if (CollectionUtils.isNotEmpty(locatedVariantList)) {
            List<Variants_80_4> variants = new ArrayList<>();
            for (LocatedVariant locatedVariant : locatedVariantList) {
                List<Variants_80_4> foundVariants = daoBean.getVariants_80_4_DAO().findByLocatedVariantId(locatedVariant.getId());
                assertTrue(CollectionUtils.isNotEmpty(foundVariants));
                variants.addAll(foundVariants);
            }
            results.addAll(binResultsFinalDiagnosticFactory.findHGMDAlmostCertainlyBenign(diagnosticBinningJob, variants));
        }
        assertTrue(CollectionUtils.isNotEmpty(results));
        logger.info("results.size(): {}", results.size());
        assertTrue(results.size() == 318);
        results.forEach(a -> logger.info(a.getLocatedVariant().getId().toString()));
    }

}