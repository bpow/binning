package org.renci.binning.diagnostic.test.ws;

import java.util.List;

import javax.ws.rs.core.Response;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.binning.core.BinningExecutorService;
import org.renci.binning.core.diagnostic.DiagnosticBinningJobInfo;
import org.renci.binning.dao.BinningDAOBeanService;
import org.renci.binning.dao.BinningDAOException;
import org.renci.binning.dao.clinbin.model.DX;
import org.renci.binning.dao.clinbin.model.DiagnosticBinningJob;
import org.renci.binning.diagnostic.test.executor.DiagnosticTestTask;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class DiagnosticTestServiceImpl implements DiagnosticTestService {

    private static final Logger logger = LoggerFactory.getLogger(DiagnosticTestServiceImpl.class);

    private BinningDAOBeanService binningDAOBeanService;

    private BinningExecutorService binningExecutorService;

    public DiagnosticTestServiceImpl() {
        super();
    }

    @Override
    public Response submit(DiagnosticBinningJobInfo info) {
        logger.debug("ENTERING submit(DiagnosticBinningJobInfo)");
        logger.info(info.toString());
        DiagnosticBinningJob binningJob = new DiagnosticBinningJob();
        try {
            binningJob.setStudy("TEST");
            binningJob.setGender(info.getGender());
            binningJob.setParticipant(info.getParticipant());
            binningJob.setListVersion(info.getListVersion());
            binningJob.setStatus(binningDAOBeanService.getDiagnosticStatusTypeDAO().findById("Requested"));
            DX dx = binningDAOBeanService.getDXDAO().findById(info.getDxId());
            binningJob.setDx(dx);
            List<DiagnosticBinningJob> foundBinningJobs = binningDAOBeanService.getDiagnosticBinningJobDAO().findByExample(binningJob);
            if (CollectionUtils.isNotEmpty(foundBinningJobs)) {
                binningJob = foundBinningJobs.get(0);
            } else {
                binningJob.setId(binningDAOBeanService.getDiagnosticBinningJobDAO().save(binningJob));
            }
            info.setId(binningJob.getId());
            logger.info(binningJob.toString());

            binningExecutorService.getExecutor().submit(new DiagnosticTestTask(binningJob.getId()));

        } catch (BinningDAOException e) {
            logger.error(e.getMessage(), e);
        }
        return Response.ok(info).build();
    }

    public BinningExecutorService getBinningExecutorService() {
        return binningExecutorService;
    }

    public void setBinningExecutorService(BinningExecutorService binningExecutorService) {
        this.binningExecutorService = binningExecutorService;
    }

    public BinningDAOBeanService getBinningDAOBeanService() {
        return binningDAOBeanService;
    }

    public void setBinningDAOBeanService(BinningDAOBeanService binningDAOBeanService) {
        this.binningDAOBeanService = binningDAOBeanService;
    }

}
