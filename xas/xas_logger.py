import logging
import logging.handlers

def add_new_print_to_logger(logger, print_func):
    info_func = logger.info

    def new_info_func(msg):
        info_func(msg)
        print_func(msg)

    logger.info = new_info_func



def get_logger(print_func=None):
    # Setup beamline specifics:
    beamline_gpfs_path = '/nsls2/xf08id'

    logger = logging.getLogger('xas_logger')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger.propagate = False

    # if print_func is not None:
    #     add_new_print_to_logger(logger, print_func)

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