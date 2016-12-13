package org.renci.binning.incidental.ncgenes.ws;

import java.util.List;

import javax.ws.rs.core.Response;

import org.apache.commons.collections4.CollectionUtils;
import org.renci.binning.core.BinningExecutorService;
import org.renci.binning.core.incidental.IncidentalBinningJobInfo;
import org.renci.binning.dao.BinningDAOBeanService;
import org.renci.binning.dao.BinningDAOException;
import org.renci.binning.dao.clinbin.model.IncidentalBinX;
import org.renci.binning.dao.clinbin.model.IncidentalBinningJob;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class IncidentalNCGenesServiceImpl implements IncidentalNCGenesService {

    private static final Logger logger = LoggerFactory.getLogger(IncidentalNCGenesServiceImpl.class);

    private BinningDAOBeanService binningDAOBeanService;

    private BinningExecutorService binningExecutorService;

    public IncidentalNCGenesServiceImpl() {
        super();
    }

    @Override
    public Response submit(IncidentalBinningJobInfo info) {
        logger.debug("ENTERING submit(IncidentalBinningJobInfo)");
        IncidentalBinningJob binningJob = new IncidentalBinningJob();
        try {
            binningJob.setStudy("GS");
            binningJob.setGender(info.getGender());
            binningJob.setParticipant(info.getParticipant());
            binningJob.setListVersion(info.getListVersion());
            binningJob.setStatus(binningDAOBeanService.getIncidentalStatusTypeDAO().findById("Requested"));
            IncidentalBinX incidentalBin = binningDAOBeanService.getIncidentalBinXDAO().findById(info.getIncidentalBinId());
            binningJob.setIncidentalBinX(incidentalBin);
            List<IncidentalBinningJob> foundBinningJobs = binningDAOBeanService.getIncidentalBinningJobDAO().findByExample(binningJob);
            if (CollectionUtils.isNotEmpty(foundBinningJobs)) {
                binningJob = foundBinningJobs.get(0);
            } else {
                binningJob.setId(binningDAOBeanService.getIncidentalBinningJobDAO().save(binningJob));
            }
            logger.info(binningJob.toString());

            // DiagnosticGSTask task = new DiagnosticGSTask();
            // task.setBinningDAOBeanService(binningDAOBeanService);
            // task.setBinningJob(binningJob);
            // binningExecutorService.getExecutor().submit(task);

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
