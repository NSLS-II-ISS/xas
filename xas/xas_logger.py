import logging
import logging.handlers


def get_logger():
    # Setup beamline specifics:
    beamline_gpfs_path = '/nsls2/xf08id'

    logger = logging.getLogger('xas_logger')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    # only add handlers if not added before
    if not len(logger.handlers):
        logger.setLevel(logging.DEBUG)
        # Write DEBUG and INFO messages to /var/log/data_processing_worker/debug.log.
        debug_file = logging.handlers.RotatingFileHandler(
            beamline_gpfs_path + '/log/data_processing_debug.log',
            maxBytes=10000000, backupCount=9)
        debug_file.setLevel(logging.DEBUG)
        debug_file.setFormatter(formatter)
        logger.addHandler(debug_file)

        # Write INFO messages to /var/log/data_processing_worker/info.log.
        info_file = logging.handlers.RotatingFileHandler(
            beamline_gpfs_path + '/log/data_processing.log',
            maxBytes=10000000, backupCount=9)
        info_file.setLevel(logging.INFO)
        info_file.setFormatter(formatter)
        logger.addHandler(info_file)


    return logger