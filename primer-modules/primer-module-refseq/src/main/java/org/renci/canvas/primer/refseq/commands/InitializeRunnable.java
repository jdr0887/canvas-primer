package org.renci.canvas.primer.refseq.commands;

import java.util.Arrays;
import java.util.List;

import org.renci.canvas.dao.CANVASDAOBeanService;
import org.renci.canvas.dao.refseq.model.GroupingType;
import org.renci.canvas.dao.refseq.model.LocationType;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class InitializeRunnable implements Runnable {

    private static final Logger logger = LoggerFactory.getLogger(PersistAction.class);

    private CANVASDAOBeanService canvasDAOBeanService;

    public InitializeRunnable(CANVASDAOBeanService canvasDAOBeanService) {
        super();
        this.canvasDAOBeanService = canvasDAOBeanService;
    }

    @Override
    public void run() {
        logger.debug("ENTERING run()");

        try {
            logger.info("Loading GroupingTypes");
            List<String> gtList = Arrays.asList("single", "order", "join");
            List<GroupingType> groupingTypeList = canvasDAOBeanService.getGroupingTypeDAO().findAll();
            for (String gt : gtList) {
                if (!groupingTypeList.stream().filter(a -> a.getId().equals(gt)).findAny().isPresent()) {
                    GroupingType groupingType = new GroupingType(gt);
                    canvasDAOBeanService.getGroupingTypeDAO().save(groupingType);
                }
            }

            logger.info("Loading LocationTypes");
            List<String> ltList = Arrays.asList("UTR", "UTR-5", "UTR-3", "intron", "exon", "intergenic", "potential RNA-editing site",
                    "intron/exon boundary");
            List<LocationType> locationTypeList = canvasDAOBeanService.getLocationTypeDAO().findAll();
            for (String lt : ltList) {
                if (!locationTypeList.stream().filter(a -> a.getId().equals(lt)).findAny().isPresent()) {
                    LocationType locationType = new LocationType(lt);
                    canvasDAOBeanService.getLocationTypeDAO().save(locationType);
                }
            }
        } catch (Exception e) {
            logger.error(e.getMessage(), e);
        }
    }

}
