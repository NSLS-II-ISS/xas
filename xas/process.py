import gc
import numpy as np
import os
import time as ttime

from .bin import bin
from .file_io import (
    load_dataset_from_files,
    create_file_header,
    validate_file_exists,
    validate_path_exists,
    save_interpolated_df_as_file,
    find_e0,
    stepscan_remove_offsets,
    stepscan_normalize_xs,
    combine_xspress3_channels,
    filter_df_by_valid_keys,
    save_primary_df_as_file,
    save_extended_data_as_file,
    dump_tiff_images,
)
from .db_io import (
    load_apb_dataset_from_db,
    translate_apb_dataset,
    load_apb_trig_dataset_from_db,
    load_xs3_dataset_from_db,
    load_pil100k_dataset_from_db,
)
from .interpolate import interpolate

from xas.file_io import _shift_root
from xas.db_io import update_header_start
from .xas_logger import get_logger

from .vonhamos import process_von_hamos_scan


def process_interpolate_bin(
    doc,
    db,
    draw_func_interp=None,
    draw_func_bin=None,
    cloud_dispatcher=None,
    print_func=None,
    dump_to_tiff=False,
):
    if "experiment" in db[doc["run_start"]].start.keys():
        uid = doc["run_start"]
        process_interpolate_bin_from_uid(
            uid,
            db,
            draw_func_interp=draw_func_interp,
            draw_func_bin=draw_func_bin,
            cloud_dispatcher=cloud_dispatcher,
            print_func=print_func,
            dump_to_tiff=dump_to_tiff,
        )


def process_interpolate_bin_from_uid(
    uid,
    db,
    draw_func_interp=None,
    draw_func_bin=None,
    cloud_dispatcher=None,
    print_func=None,
    dump_to_tiff=False,
    save_interpolated_file=True,
    update_start=None,
):
    logger = get_logger()
    (
        hdr,
        primary_df,
        extended_data,
        comments,
        path_to_file,
        file_list,
        data_kind,
    ) = get_processed_df_from_uid(
        uid,
        db,
        logger=logger,
        draw_func_interp=draw_func_interp,
        draw_func_bin=draw_func_bin,
        print_func=print_func,
        save_interpolated_file=save_interpolated_file,
        update_start=update_start,
    )

    save_primary_df_as_file(path_to_file, primary_df, comments)

    try:
        save_extended_data_as_file(path_to_file, extended_data, data_kind=data_kind)
    except Exception as e:
        print(e)
        pass

    if dump_to_tiff:
        if extended_data is not None:
            tiff_files = dump_tiff_images(path_to_file, primary_df, extended_data)
            print(f" >>>>>>>>>>Tiff file {tiff_files}")
            file_list += tiff_files

    # dfgd
    try:
        if cloud_dispatcher is not None:
            year, cycle, proposal = (
                hdr.start["year"],
                hdr.start["cycle"],
                hdr.start["PROPOSAL"],
            )
            for f in file_list:
                cloud_dispatcher.load_to_dropbox(
                    f,
                    year=year,
                    cycle=cycle,
                    proposal=proposal,
                )
                logger.info(
                    f"({ttime.ctime()}) Sending data to the cloud successful for {path_to_file}"
                )
    except Exception as e:
        logger.info(
            f"({ttime.ctime()}) Sending data to the cloud failed for {path_to_file}"
        )
        raise e

    clear_db_cache(db)


_legacy_experiment_reg = {"fly_energy_scan_pil100k": "fly_scan"}


def get_processed_df_from_uid(
    uid,
    db,
    logger=None,
    draw_func_interp=None,
    draw_func_bin=None,
    print_func=None,
    save_interpolated_file=True,
    update_start=None,
):
    if print_func is None:
        print_func = print
    if logger is None:
        logger = get_logger()

    hdr = db[uid]
    if update_start is not None:
        hdr = update_header_start(hdr, update_start)
    experiment = hdr.start["experiment"]
    if experiment in _legacy_experiment_reg.keys():
        experiment = _legacy_experiment_reg[experiment]
    comments = create_file_header(hdr)
    path_to_file = hdr.start["interp_filename"]
    path_to_file = _shift_root(path_to_file)
    validate_path_exists(path_to_file)
    path_to_file = validate_file_exists(path_to_file, file_type="interp")
    e0 = find_e0(hdr)
    data_kind = "default"
    file_list = []

    if experiment == "fly_scan":
        # path_to_file = validate_file_exists(path_to_file, file_type='interp')
        stream_names = hdr.stream_names
        try:
            # default detectors
            (
                apb_df,
                energy_df,
                energy_offset,
            ) = load_apb_dataset_from_db(db, uid)
            raw_dict = translate_apb_dataset(apb_df, energy_df, energy_offset)

            for stream_name in stream_names:
                if stream_name == "pil100k_stream":
                    apb_trigger_pil100k_timestamps = load_apb_trig_dataset_from_db(
                        db,
                        uid,
                        use_fall=True,
                        stream_name="apb_trigger_pil100k",
                    )
                    pil100k_dict = load_pil100k_dataset_from_db(
                        db,
                        uid,
                        apb_trigger_pil100k_timestamps,
                    )
                    raw_dict = {**raw_dict, **pil100k_dict}

                elif stream_name == "xs_stream":
                    apb_trigger_xs_timestamps = load_apb_trig_dataset_from_db(
                        db,
                        uid,
                        stream_name="apb_trigger_xs",
                    )
                    xs3_dict = load_xs3_dataset_from_db(
                        db, uid, apb_trigger_xs_timestamps
                    )
                    raw_dict = {**raw_dict, **xs3_dict}

            logger.info(
                f"({ttime.ctime()}) Loading file successful for UID {uid}/{path_to_file}"
            )
        except Exception as e:
            logger.info(
                f"({ttime.ctime()}) Loading file failed for UID {uid}/{path_to_file}"
            )
            raise e
        try:
            interpolated_df = interpolate(raw_dict)
            logger.info(
                f"({ttime.ctime()}) Interpolation successful for {path_to_file}"
            )
            if save_interpolated_file:
                save_interpolated_df_as_file(path_to_file, interpolated_df, comments)
        except Exception as e:
            logger.info(f"({ttime.ctime()}) Interpolation failed for {path_to_file}")
            raise e

        try:
            if e0 > 0:
                processed_df = bin(interpolated_df, e0)
                (path, extension) = os.path.splitext(path_to_file)
                path_to_file = path + ".dat"
                logger.info(f"({ttime.ctime()}) Binning successful for {path_to_file}")

                if draw_func_interp is not None:
                    draw_func_interp(interpolated_df)
                if draw_func_bin is not None:
                    draw_func_bin(processed_df, path_to_file)
            else:
                print(f"({ttime.ctime()}) Energy E0 is not defined")
        except Exception as e:
            logger.info(f"({ttime.ctime()}) Binning failed for {path_to_file}")
            raise e

        # save_binned_df_as_file(path_to_file, processed_df, comments)

    elif (experiment == "step_scan") or (experiment == "collect_n_exposures"):
        # path_to_file = validate_file_exists(path_to_file, file_type='interp')
        df = stepscan_remove_offsets(hdr)
        df = stepscan_normalize_xs(df)
        processed_df = filter_df_by_valid_keys(df)
        # df_processed = combine_xspress3_channels(df)

    else:
        return

    processed_df = combine_xspress3_channels(processed_df)

    (
        primary_df,
        extended_data,
    ) = split_df_data_into_primary_and_extended(processed_df)

    ### WIP
    if "spectrometer" in hdr.start.keys():
        if hdr.start["spectrometer"] == "von_hamos":
            (
                extended_data,
                comments,
                file_paths,
            ) = process_von_hamos_scan(
                primary_df,
                extended_data,
                comments,
                hdr,
                path_to_file,
                db=db,
            )
            data_kind = "von_hamos"
            file_list = file_paths
            # save_vh_scan_to_file(path_to_file, vh_scan, comments)
    file_list.append(path_to_file)
    return (
        hdr,
        primary_df,
        extended_data,
        comments,
        path_to_file,
        file_list,
        data_kind,
    )


def get_df_and_metadata_from_db(
    uid,
    db,
    logger=None,
    draw_func_interp=None,
    draw_func_bin=None,
    print_func=None,
    save_interpolated_file=True,
    update_start=None,
):
    """Summary
    
    Parameters
    ----------
    uid : TYPE
        Description
    db : TYPE
        Description
    logger : None, optional
        Description
    draw_func_interp : None, optional
        Description
    draw_func_bin : None, optional
        Description
    print_func : None, optional
        Description
    save_interpolated_file : bool, optional
        Description
    update_start : None, optional
        Description

    Returns
    -------
    dict
        A dictionary with two keys: "data" and "metadata". The data is simply a
        Pandas dataframe and the metadata is a python dictionary.
    """

    # TODO
    raise NotImplementedError


def split_df_data_into_primary_and_extended(df_orig):
    sec_cols = []
    for c in df_orig.columns:
        if len(df_orig[c][0].shape) > 0:
            sec_cols.append(c)

    if len(sec_cols) > 0:
        extended_data = {}
        df = df_orig.copy()
        for col in sec_cols:
            v = df.pop(col)
            extended_data[col] = np.array([i for i in v])

    else:
        df = df_orig
        extended_data = None
    return df, extended_data


def clear_db_cache(db):
    db._catalog._entries.cache_clear()
    gc.collect()


def process_interpolate_only(doc, db):
    if "experiment" in db[doc["run_start"]].start.keys():
        if db[doc["run_start"]].start["experiment"] == "fly_energy_scan":
            raw_df = load_dataset_from_files(db, doc["run_start"])
            interpolated_df = interpolate(raw_df)
            return interpolated_df


def process_interpolate_unsorted(uid, db):
    raw_df = load_dataset_from_files(db, uid)
    interpolated_df = interpolate(raw_df, sort=False)
    return interpolated_df
