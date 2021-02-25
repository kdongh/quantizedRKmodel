import matplotlib.pyplot as plt
import numpy as np

import multi_function


class visualize:

    def __init__(self, si_arguments, arguments, parameter):
        self.number_of_sites = arguments['number_of_sites']

        self.time, self.pulse, self.data = multi_function.td_schrodinger(si_arguments, arguments)
        self.make_plot(parameter)

    def data_divide_site(self):
        result = []
        self.data = np.array(self.data)
        for i in range(self.number_of_sites + 1):
            result.append(self.data[..., i])

        return result

    def make_plot(self, title):
        data = self.data_divide_site()
        for i in range(self.number_of_sites + 1):
            plt.plot(self.time, data[i], label=i)

        plt.plot(self.time, self.pulse, label='pulse')

        plt.legend()
        plt.title(title)
        plt.ylim([-0.2,1.2])
        plt.savefig(str(title) + '.png')
        #plt.show()