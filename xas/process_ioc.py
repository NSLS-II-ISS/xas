from textwrap import dedent

from caproto import ChannelType
from caproto._dbr import _LongStringChannelType
from caproto.server import PVGroup, ioc_arg_parser, pvproperty, run
import time as ttime
import numpy as np
from ophyd import Component as Cpt, Device, EpicsSignal, Kind


print('#####################################################################')
print('###                   LAUNCHING PROCESSING IOC                    ###')
print('#####################################################################')

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
#
#     @pv_sp.startup
#     async def pv_sp(self, instance, async_lib):
#         """This is a startup hook which periodically updates the value."""
#
#         while True:
#             await async_lib.sleep(self.dwell)
#             if self.pv_data_available == 1:
#
#
#
#             if self.go.value == 1:
#                 if self.time_start is None:
#                     self.time_start = ttime.time()
#
#                 if self.pause.value == 0:
#
#                     await self.time_elapsed.write(ttime.time() - self.time_start)
#                     dt = self.time_elapsed.value - self.time_paused.value
#
#                     if dt < np.max(self.tprog.value):
#                         pv_sp = np.interp(dt, self.tprog.value, self.pvprog.value)
#                     else:
#                         pv_sp = self.pvprog.value[-1]
#
#                     if self._previous_elapsed_time is not None:
#                         rate = (pv_sp - self.pv_sp.value) / (dt - self._previous_elapsed_time) * 60
#
#                     self._previous_elapsed_time = dt
#                     await instance.write(value=pv_sp)
#                 else:
#                     rate = 0
#
#             else:
#                 self.time_start = None
#                 await self.time_elapsed.write(0.0)
#                 await self.time_paused.write(0.0)
#                 rate = 0
#
#
#
#     @pid_enable_name.startup
#     async def pid_enable_name(self, instance, async_lib):
#         while self.pid_enable_name.value == '':
#             await async_lib.sleep(1)
#
#         self.pid_enable = EpicsSignal(self.pid_enable_name.value, name='pid_enable')
#         def subscription(value, **kwargs):
#             if value == 0:
#                 if self.pid_output is not None:
#                     if self.pid_output.get() != 0:
#                         self.pid_output.put(0)
#         self.pid_enable.subscribe(subscription)
#
#
#     @pid_output_name.startup
#     async def pid_output_name(self, instance, async_lib):
#         while (self.pid_output_name.value == ''):
#             await async_lib.sleep(1)
#
#         self.pid_output = EpicsSignal(self.pid_output_name.value, name='pid_output')
#         def subscription(value, **kwargs):
#             if value != 0:
#                 if self.pid_enable is not None:
#                     if self.pid_enable.get() == 0:
#                         self.pid_output.put(0)
#         self.pid_output.subscribe(subscription)
#
#
#     @safety_timer.startup
#     async def safety_timer(self, instance, async_lib):
#         while self.safety_timer.value < self.safety_thresh.value:
#             safety_timer = self.safety_timer.value + 1
#             await instance.write(value=safety_timer)
#             await async_lib.sleep(1)
#         if self.pid_enable is not None:
#             self.pid_enable.put(0)
#
#
#     @time_paused.startup
#     async def time_paused(self, instance, async_lib):
#         while True:
#             cur_time = ttime.time()
#             if self.pause.value == 1:
#                 while True:
#                     if self.pause.value == 0:
#                         break
#                     await async_lib.sleep(self.dwell.value/2)
#                 time_paused = self.time_paused.value + (ttime.time() - cur_time)
#                 await instance.write(value=time_paused)
#
#             await async_lib.sleep(self.dwell.value)
#
#
#
#
#
#
#
#
# if __name__ == '__main__':
#     ioc_options, run_options = ioc_arg_parser(
#         default_prefix='XF:08IDB-Ramping:',
#         desc=dedent(ProcessingIOC.__doc__))
#     ioc = ProcessingIOC(**ioc_options)
#     run(ioc.pvdb, **run_options)