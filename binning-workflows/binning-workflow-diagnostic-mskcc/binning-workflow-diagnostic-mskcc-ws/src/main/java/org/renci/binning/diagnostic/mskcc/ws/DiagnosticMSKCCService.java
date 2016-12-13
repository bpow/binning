package org.renci.binning.diagnostic.mskcc.ws;

import javax.ws.rs.Consumes;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import org.renci.binning.core.diagnostic.DiagnosticBinningJobInfo;

@Path("/DiagnosticMSKCCService/")
@Consumes(MediaType.APPLICATION_JSON)
@Produces(MediaType.APPLICATION_JSON)
public interface DiagnosticMSKCCService {

    @POST
    @Path("/submit")
    public Response submit(DiagnosticBinningJobInfo info);

}
