package org.renci.canvas.primer.gr.commands;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "genome-reference", name = "persist", description = "Persist")
@Service
public class PersistAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(PersistAction.class);

    @Reference
    private CANVASDAOBeanService canvasDAOBeanService;

    public PersistAction() {
        super();
    }

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");

        Executors.newSingleThreadExecutor().submit(new PersistRunnable(canvasDAOBeanService));
        return null;
    }
}
