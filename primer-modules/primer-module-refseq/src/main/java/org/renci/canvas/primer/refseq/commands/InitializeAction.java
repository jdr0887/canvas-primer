package org.renci.canvas.primer.refseq.commands;

import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.renci.canvas.dao.CANVASDAOBeanService;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

@Command(scope = "refseq", name = "init", description = "Initialize")
@Service
public class InitializeAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(PersistAction.class);

    @Reference
    private CANVASDAOBeanService canvasDAOBeanService;

    public InitializeAction() {
        super();
    }

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");
        long start = System.currentTimeMillis();
        Executors.newSingleThreadExecutor().execute(new InitializeRunnable(canvasDAOBeanService));
        long end = System.currentTimeMillis();
        logger.info("duration: {} seconds", (end - start) / 1000D);
        return null;

    }

}
