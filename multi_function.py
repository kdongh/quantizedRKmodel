import multiprocessing as multi
import numpy as np
import csv

import runge_kutta


def number_to_combination(number, number_of_sites):
    text = "{0:b}".format(number).zfill(number_of_sites)
    result = []
    for i in range(number_of_sites):
        result.append(int(text[i]))
    return result


def combination_to_number(combination):
    str_combination = ''.join(list(map(str, combination)))
    result = int(str_combination, base=2)
    return result


def wave_to_population(wavefunction):
    result = []

    for i in range(len(wavefunction)):
        result.append(np.abs(wavefunction[i]))

    return result


def wave_to_real(wavefunction):
    result = []

    for i in range(len(wavefunction)):
        result.append(np.real(wavefunction[i]))

    return result


def generate_basis(number_of_sites):
    '''
    generate basis for Hamiltonian
    number_of_sites + 1 = index 0 means no particle state and sites are counting from 1 not 0
    :param number_of_sites: int
    :return: numpy array
    '''
    basis = np.zeros(number_of_sites + 1, dtype=complex)
    return basis


def particle_energy(start_int, end_int, wavefunction_in, wavefunction_out, particle_frequency):
    '''
    particle energy return hbar = 1 case
    :param start_int: int
    :param end_int: int
    :param wavefunction_in: numpy array
    :param wavefunction_out: numpy array
    :param particle_frequency: float
    :return: float
    '''
    caliblation = 1 / 2
    for i_site in range(start_int, end_int):
        if i_site:
            wavefunction_out[i_site] += particle_frequency * wavefunction_in[i_site]
            # print(i_site,wavefunction_out[i_site])
        wavefunction_out[i_site] += particle_frequency * caliblation * wavefunction_in[i_site]


def hopping(start_int, end_int, number_of_sites, wavefunction_in, wavefunction_out, hopping_t):
    '''
    hopping energy return to wavefunction_out with hbar = 1
    :param start_int: int
    :param end_int: int
    :param number_of_sites: int
    :param wavefunction_in: numpy array
    :param wavefunction_out: numpy array
    :param particle_frequency: float
    :return:
    '''
    hopping_energy = -hopping_t
    for i_site in range(start_int, end_int):
        if i_site != 0:
            if i_site == 1:  # left boundary hopping
                energy = hopping_energy * wavefunction_in[1]
                wavefunction_out[number_of_sites] += energy  # to left
                wavefunction_out[2] += energy  # to right
            elif i_site == number_of_sites:  # right boundary hopping
                energy = hopping_energy * wavefunction_in[number_of_sites]
                wavefunction_out[number_of_sites - 1] += energy  # to left
                wavefunction_out[1] += energy  # to right
            else:
                energy = hopping_energy * wavefunction_in[i_site]
                wavefunction_out[i_site - 1] += energy
                wavefunction_out[i_site + 1] += energy


def pulse_amp(time, arguments, step_pulse):
    pulse_amp = arguments['pulse_amp']
    pulse_frequency = arguments['pulse_frequency']
    pulse_delay = arguments['pulse_delay_time']
    pulse_average_time = arguments['pulse_average_time']
    pulse_length = arguments['pulse_length']

    delayed_time = time - pulse_delay
    pulse_start = pulse_delay - pulse_length/2.0
    pulse_end = pulse_delay + pulse_length/2.0
    if step_pulse:
        if time < pulse_end or time > pulse_start:
            result = pulse_amp * np.exp(1j * pulse_frequency * delayed_time) * np.exp(
            -(delayed_time ** 2) / (pulse_average_time ** 2))
        else:
            result = 0
    else:
        result = pulse_amp * np.exp(1j * pulse_frequency * delayed_time) * np.exp(
            -(delayed_time ** 2) / (pulse_average_time ** 2))

    return result


def pulse(time, wavefunction_in, wavefunction_out, arguments, step_pulse):
    amp = pulse_amp(time, arguments, step_pulse)

    wavefunction_out[1] += amp * wavefunction_in[0]
    wavefunction_out[0] += np.conjugate(amp) * wavefunction_in[1]


def imaginary(wavefunction):
    return -1j * wavefunction


def normalization(wavefunction):
    buff = 0
    for i in range(len(wavefunction)):
        buff += np.abs(wavefunction[i]) ** 2
    buff = np.sqrt(buff)
    if not buff:
        buff = 1
    return wavefunction / buff


def checking_norm(wavefunction):
    buff = 0
    for i in range(len(wavefunction)):
        buff += np.abs(wavefunction[i]) ** 2

    return buff


def total_hamiltonian(time, wavefunction_in, arguments, step_pulse):
    number_of_sites = arguments['number_of_sites']
    particle_frequency = arguments['particle_frequency']
    hopping_t = arguments['hopping_t']
    buffer_wavefunction = np.zeros(number_of_sites + 1, dtype=complex)
    particle_energy(0, number_of_sites + 1, wavefunction_in, buffer_wavefunction, particle_frequency)
    # print('after energy = ',buffer_wavefunction)
    hopping(0, number_of_sites + 1, number_of_sites, wavefunction_in, buffer_wavefunction, hopping_t)
    pulse(time, wavefunction_in, buffer_wavefunction, arguments, step_pulse)
    # print('after pulse = ',buffer_wavefunction)

    # buffer_wavefunction = normalization(buffer_wavefunction)
    # print('after norm = ',buffer_wavefunction)

    return buffer_wavefunction


def td_schrodinger(si_arguments, arguments, step_pulse=False, outfile_name='none.csv', csv_out=False):
    si_dt = si_arguments['time_gap']

    record_integer = 0.5 / si_dt

    end_time = arguments['time_end']
    dt = arguments['time_gap']

    number_of_sites = arguments['number_of_sites']

    total_number_of_time_integer = int(end_time // dt)

    wavefunction = generate_basis(number_of_sites)
    wavefunction[0] = 1  # for empty state
    wavefunction[1] = 0

    time_list = []
    pulse_list = []
    result_list = []
    norm_list = []
    log = False

    if csv_out:
        out_file_name = outfile_name + '.csv'
        out_file = open(out_file_name, 'w', newline='')
        out_writer = csv.writer(out_file)

    for i_time in range(total_number_of_time_integer):
        time = i_time * dt
        # print('before runge = ',wavefunction)
        wavefunction = runge_kutta.runge_kutta(time, wavefunction, total_hamiltonian, arguments, log, step_pulse)
        # print('before Norm = ',wavefunction)
        wavefunction = normalization(wavefunction)
        # print('after Norm = ',wavefunction)
        if i_time % record_integer == 0:
            print(i_time, ' / ', total_number_of_time_integer)
            time_list.append(i_time * si_dt)
            pulse_result = np.real(pulse_amp(time, arguments, step_pulse))
            pulse_list.append(pulse_result)
            # population_result = wave_to_population(wavefunction)
            real_result = wave_to_real(wavefunction)
            # result_list.append(population_result)
            norm_list.append(checking_norm(wavefunction))
            if csv_out:
                real_result.insert(0, pulse_result)
                real_result.insert(0, i_time * si_dt)
                out_writer.writerow(real_result)

    if csv_out:
        out_file.close()
    return time_list, pulse_list, result_list, norm_list


def function_to_multi(time, wavefunction_in, wavefunction_out, arguments):
    number_of_sites = int(arguments['number_of_sites'])
    number_of_processes = int(arguments['number_of_processes'])
    particle_frequency = int(arguments['particle_frequency'])
    total_number_of_case = 2 ** number_of_sites

    case_per_process = total_number_of_case // number_of_processes
    rest_case = total_number_of_case - 1

    for i_multi in range(number_of_processes):
        start_point = i_multi * case_per_process
        end_point = start_point + case_per_process

        multi.Process(target=total_hamiltonian, args=(
            start_point, end_point, number_of_sites, wavefunction_in, wavefunction_out, particle_frequency))

    return 0
