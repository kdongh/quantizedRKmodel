import numpy as np
def argv_to_parameter(argv, parameter):
    parameter['number_of_processes'] = int(argv[1])

    parameter['number_of_sites'] = int(argv[2])

    parameter['time_end'] = float(argv[3])# fs
    parameter['time_gap'] = float(argv[4])# fs
    ''' frequency approximation
    parameter['particle_mass'] = argv[5]# electron mass
    parameter['rod_length'] = argv[6]# angstrom
    parameter['torsion_constant'] = argv[7]# meV

    parameter['g_acceleration'] = 9.764 # m/s^2
    '''
    parameter['particle_frequency'] = float(argv[5])# THz
    parameter['hopping_t'] = float(argv[6])# THz

    parameter['pulse_amp'] = float(argv[7])
    parameter['pulse_frequency'] = float(argv[8])# THz
    parameter['pulse_delay_time'] = float(argv[9])# fs
    parameter['pulse_average_time'] = float(argv[10])# fs
    parameter['pulse_length'] = float(argv[11])# fs


class si_to_au:

    def __init__(self, si, au):

        au['number_of_processes'] = si['number_of_processes']
        au['number_of_sites'] = si['number_of_sites']

        au['time_end'] = self.time_fs(si['time_end'])
        au['time_gap'] = self.time_fs(si['time_gap'])

        #au['particle_mass'] = self.mass_electron(si['particle_mass'])
        #au['rod_length'] = self.length_angstrom(si['rod_length'])
        #au['torsion_constant'] = self.energy_meV(si['torsion_constant'])

        au['particle_frequency'] = self.frequency(si['particle_frequency'])
        au['hopping_t'] = self.frequency(si['hopping_t'])

        au['pulse_amp'] = si['pulse_amp']
        au['pulse_frequency'] = self.frequency(si['pulse_frequency'])
        au['pulse_delay_time'] = self.time_fs(si['pulse_delay_time'])
        au['pulse_average_time'] = self.time_fs(si['pulse_average_time'])
        au['pulse_length'] = self.time_fs(si['pulse_length'])

    def energy_meV(self, meV):# 1au = 27.211386245 eV
        au = meV / 27211.386245
        return au

    def time_fs(self, fs):# 1au = 0.02418884326 fs
        au = fs / 0.02418884326
        return au

    def mass_electron(self, ele):# 1au = 9.1093837015 * 10^-31kg
        au = ele
        return au

    def length_angstrom(self, a):# 1au = 0.529177210903 angstrom
        au = a / 0.529177210903
        return au

    def frequency(self, THz):# 1THz = 0.001 / fs = 0.001 * 0.02418884326 / au
        freq = 2 * np.pi * THz * 0.001 * 0.0241888436
        return freq
