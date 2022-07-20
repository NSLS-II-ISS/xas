from textwrap import dedent

from caproto import ChannelType
from xas.process import process_interpolate_bin_from_uid
from xas.handlers import get_iss_db
from caproto.server import PVGroup, ioc_arg_parser, pvproperty, run

import numpy as np
from ophyd import Component as Cpt, Device, EpicsSignal, Kind

print('#####################################################################')
print('###                   LAUNCHING PROCESSING IOC                    ###')
print('#####################################################################')

db = get_iss_db()

class ProcessingIOC(PVGroup):
    """
    An IOC for executing ramps of pv_setpoints

    """
    pv_data_available = pvproperty(
        value=0,
        doc='pv data available',
    )

    pv_uid = pvproperty(
        value='',
        dtype=ChannelType.STRING,
        doc='pv uid',
    )

    dwell = 0.5
    uid_list = []

    @pv_uid.startup
    async def pv_uid(self, instance, async_lib):
        """This is a startup hook which periodically updates the value."""
        while True:
            await async_lib.sleep(self.dwell)
            # print(f'{instance.value=}')
            if instance.value != '':
                self.uid_list.append(instance.value)
                await instance.write(value='')
                await self.pv_data_available.write(1)

    @pv_data_available.startup
    async def pv_data_available(self, instance, async_lib):
        """This is a startup hook which periodically updates the value."""
        #
        while True:
            await async_lib.sleep(self.dwell)
            if instance.value == 1:
                process_uid = self.uid_list.pop(0)

                print(f'PROCESS IOC: File received {process_uid}')
                process_interpolate_bin_from_uid(process_uid, db)

                await instance.write(0)

if __name__ == '__main__':
    ioc_options, run_options = ioc_arg_parser(
        default_prefix='XF:08IDB-Processing:',
        desc=dedent(ProcessingIOC.__doc__))
    ioc = ProcessingIOC(**ioc_options)
    run(ioc.pvdb, **run_options)