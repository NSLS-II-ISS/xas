from iss_workflows.metadata import metadata_key_dict as metadata_dict
from iss_workflows.metadata import metadata_key_match, ghs_selected_gas_key_match


def generate_file_header_from_hdr(hdr):
    output = ''
    for key, hr_key in metadata_key_match.items():
        if key in hdr.start.keys():
            value = hdr.start[key]
        elif key == 'stop_time':
            value = hdr.stop['time']
        else:
            value = 'None'
        if (key == 'time') or (key == 'stop_time'):
            value = datetime.fromtimestamp(value).strftime('%m/%d/%Y  %H:%M:%S.%f')
        output += f'# {hr_key}: {value}\n'
    return output




