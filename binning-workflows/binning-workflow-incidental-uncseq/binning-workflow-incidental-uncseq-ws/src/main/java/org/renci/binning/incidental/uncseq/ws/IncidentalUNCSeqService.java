package org.renci.binning.incidental.uncseq.ws;

import javax.ws.rs.Consumes;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import org.renci.binning.core.incidental.IncidentalBinningJobInfo;

@Path("/IncidentalUNCSeqService/")
@Consumes(MediaType.APPLICATION_JSON)
@Produces(MediaType.APPLICATION_JSON)
public interface IncidentalUNCSeqService {

    @POST
    @Path("/submit")
    public Response submit(IncidentalBinningJobInfo info);

}
