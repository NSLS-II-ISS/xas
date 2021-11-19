

def show_detector_time_traces(db, uid, fig=None):
    hdr = db[uid]
    for stream_name in hdr.stream_names:
        t = hdr.table(stream_name=stream_name, fill=True)





