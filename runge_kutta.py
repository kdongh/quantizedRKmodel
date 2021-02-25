import numpy as np


def runge_kutta(time, initial_value, function, arguments, log):
    h = arguments['time_gap']
    k1 = -1j * h * function(time, initial_value, arguments)
    k2 = -1j * h * function(time + 0.5 * h, initial_value + 0.5 * k1, arguments)
    k3 = -1j * h * function(time + 0.5 * h, initial_value + 0.5 * k2, arguments)
    k4 = -1j * h * function(time + h, initial_value + k3, arguments)

    result = initial_value + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    if log:
        print('energy = ',arguments['particle_frequency'])
        print('initial = ',initial_value)
        print('k1 = ',k1)
        print('k2 = ',k2)
        print('k3 = ',k3)
        print('k4 = ',k4)
        print('result = ',result)
    return result
