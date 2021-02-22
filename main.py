import sys

import multi_function
import parameter_setting

import visualize

parameter_si = {}
parameter_au = {}


if __name__ == '__main__':
    if len(sys.argv) != 10:
        print("Less argument error")
        sys.exit()

    parameter_setting.argv_to_parameter(sys.argv, parameter_si)
    '''
    parameter_si['time_end'] = 10
    parameter_si['time_gap'] = 0.0001

    parameter_si['particle_frequency'] = 2

    parameter_si['pulse_amp'] = 1
    parameter_si['pulse_frequency'] = 2
    parameter_si['pulse_delay_time'] = 400
    parameter_si['pulse_average_time'] = 50
    '''

    print('Start ',parameter_si)

    parameter_setting.si_to_au(parameter_si,parameter_au)
    visualize.visualize(parameter_si,parameter_au,parameter_si['pulse_average_time'])

    print('End ', parameter_si)