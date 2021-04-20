import sys

import multi_function
import parameter_setting

#import visualize

parameter_si = {}
parameter_au = {}


if __name__ == '__main__':
    if len(sys.argv) != 12:
        print("Less argument error")
        sys.exit()

    parameter_setting.argv_to_parameter(sys.argv, parameter_si)
    step_pulse_bool = False
    if parameter_si['pulse_average_time'] > 100000:
        step_pulse_bool = True
    if False:
        parameter_si['number_of_processes'] = 1
        parameter_si['number_of_sites'] = 20

        parameter_si['time_end'] = 10
        parameter_si['time_gap'] = 0.005

        parameter_si['particle_frequency'] = 3
        parameter_si['hopping_t'] = 3

        parameter_si['pulse_amp'] = 0.00001
        parameter_si['pulse_frequency'] = 1
        parameter_si['pulse_delay_time'] = 500
        parameter_si['pulse_average_time'] = 200


    print('Start ', parameter_si)

    plot_name = str(parameter_si['pulse_average_time'])+'_'+str(parameter_si['pulse_amp']) + '_' + str(parameter_si['pulse_frequency'])
    parameter_setting.si_to_au(parameter_si, parameter_au)
    #visualize.visualize(parameter_si, parameter_au, plot_name)
    multi_function.td_schrodinger(parameter_si,parameter_au,step_pulse=step_pulse_bool,outfile_name=plot_name,csv_out=True)


    print('End ', parameter_si)