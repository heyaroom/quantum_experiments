import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .randomized_benchmarking import RandomizedBenchmarking
from ...objects import Report, Job, JobTable

def exp_decay(x,a,b,p):
    y = a*p**x + b
    return y

def double_exp_decay(x,p1,p2):
    y = 1/3.*p1**x + 2/3.*p2**x
    return y

class AdjointRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.standard_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, initial_inverse=False, interleaved=interleaved)
        self.inversed_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, initial_inverse=True, interleaved=interleaved)
        self.number_of_qubit = self.standard_rb.number_of_qubit
        self.length_list = self.standard_rb.length_list

    def execute(self, take_data):
        self.standard_rb.execute(take_data)
        self.inversed_rb.execute(take_data)

    def analyze(self):
        self.standard_rb.analyze()
        self.inversed_rb.analyze()

        self.variance_list = []
        for key in self.length_list:
            variance = np.mean((np.array(self.standard_rb.data_table[key]) - np.array(self.inversed_rb.data_table[key]))**2)
            self.variance_list.append(variance)

        est_p0 = self.standard_rb.p
        est_p1 = self.standard_rb.p

        popt, pcov = curve_fit(double_exp_decay,self.length_list,self.variance_list,p0=[est_p0, est_p1])

        self.p0 = popt[0]
        self.p1 = popt[1]

        self.report = Report(name="adjoint_randomized_benchmarking")
        self.report.add_information("variance", self.variance_list)
        self.report.add_information("p0", self.p0)
        self.report.add_information("p1", self.p1)
        self.report.add_information("standard_report", self.standard_rb.report.dictionary)
        self.report.add_information("inversed_report", self.inversed_rb.report.dictionary)

    def visualize(self):
        xfit = np.linspace(0,self.length_list[-1],1001)

        plt.figure(figsize=(10,5))
        plt.subplot(121)
        srb_yfit = exp_decay(xfit,self.standard_rb.a,self.standard_rb.b,self.standard_rb.p)
        irb_yfit = exp_decay(xfit,self.inversed_rb.a,self.inversed_rb.b,self.inversed_rb.p)
        plt.plot(xfit,srb_yfit,'r-', label="Standard")
        plt.plot(xfit,irb_yfit,'b-', label="Inversed")
        plt.errorbar(x=self.length_list, y=self.standard_rb.pauli_ave, yerr=self.standard_rb.pauli_std, fmt='k.')
        plt.errorbar(x=self.length_list, y=self.inversed_rb.pauli_ave, yerr=self.inversed_rb.pauli_std, fmt='k.')
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.legend()

        plt.subplot(122)
        var_yfit = double_exp_decay(xfit, self.p0, self.p1)
        plt.plot(xfit, var_yfit, "r-")
        plt.plot(self.length_list, self.variance_list, 'k.-')
        plt.xlabel('Sequence length')
        plt.ylabel("Variance")
        plt.ylim(-0.1,1.1)
        plt.tight_layout()
        plt.show()