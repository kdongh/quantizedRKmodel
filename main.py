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
    if False:
        parameter_si['number_of_sites'] = 20

        parameter_si['time_end'] = 5000
        parameter_si['time_gap'] = 0.005

        parameter_si['particle_frequency'] = 3
        parameter_si['hopping_t'] = 3

        parameter_si['pulse_amp'] = 0.00001
        parameter_si['pulse_frequency'] = 1
        parameter_si['pulse_delay_time'] = 500
        parameter_si['pulse_average_time'] = 200

        for i in range(3):
            parameter_si['pulse_frequency'] = 2.8 + i * 0.2
            print('Start ',parameter_si)

            plot_name = str(parameter_si['pulse_amp'])+', '+str(parameter_si['pulse_frequency'])
            parameter_setting.si_to_au(parameter_si,parameter_au)
            visualize.visualize(parameter_si,parameter_au,plot_name)

            print('End ', parameter_si)

    print('Start ', parameter_si)

    plot_name = str(parameter_si['pulse_amp']) + '_' + str(parameter_si['pulse_frequency'])
    parameter_setting.si_to_au(parameter_si, parameter_au)
    visualize.visualize(parameter_si, parameter_au, plot_name)

    print('End ', parameter_si)