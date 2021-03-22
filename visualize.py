import matplotlib.pyplot as plt
import numpy as np

import multi_function


class visualize:

    def __init__(self, si_arguments, arguments, parameter):
        self.number_of_sites = arguments['number_of_sites']

        self.time, self.pulse, self.data, self.norm = multi_function.td_schrodinger(si_arguments, arguments)
        #self.make_plot(parameter)
        self.make_color_plot(parameter)

    def data_divide_site(self):
        result = []
        self.data = np.array(self.data)
        for i in range(self.number_of_sites + 1):
            result.append(self.data[..., i])

        return result

    def make_plot(self, title):
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        data = self.data_divide_site()
        for i in range(1, self.number_of_sites + 1):
            ax1.plot(self.time, data[i], label=i)
        axis = np.ones(len(self.time))
        ax2.plot(self.time, self.pulse, label='pulse')
        #ax1.plot(self.time,axis,lw=1,linestyle='dashed')
        #ax1.plot(self.time,self.norm,lw=1,linestyle='dashed',label='norm')

        ax1.legend()
        plt.title(title)
        #ax1.set_ylim([-0.5,1.5])
        ax2.set_ylim([-0.001,0.001])
        #plt.savefig(str(title) + '.png')
        plt.show()

    def make_color_plot(self, title):
        self.data = np.array(self.data)[0:,1:]
        shape = self.data.shape
        rotated_data = np.zeros((shape[1], shape[0]))
        for i in range(len(self.data)):
            for j in range(len(self.data[0])):
                rotated_data[j][i] = self.data[i][j]
        plt.pcolormesh(rotated_data)
        plt.colorbar()
        plt.title(title)
        plt.savefig(str(title) + '.png')
        #plt.show()