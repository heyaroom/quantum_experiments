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

class InterleavedRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.standard_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, interleaved=None)
        self.interleaved_rb = RandomizedBenchmarking(circuit, qubit_index, group, sequence_list, seed, interleaved)

    def execute(self, take_data):
        self.standard_rb.execute(take_data)
        self.interleaved_rb.execute(take_data)

    def analyze(self):
        self.standard_rb.analyze()
        self.interleaved_rb.analyze()

        number_of_qubit = self.standard_rb.number_of_qubit
        new_p = self.interleaved_rb.p/self.standard_rb.p
        self.fidelity = (1 + (2**number_of_qubit - 1)*new_p)/(2**number_of_qubit)

        self.report = Report(name="interleaved_randomized_benchmarking")
        self.report.add_information("average gate fidelty", self.fidelity)
        self.report.add_information("standard_report", self.standard_rb.report.dictionary)
        self.report.add_information("interleaved_report", self.interleaved_rb.report.dictionary)     

    def visualize(self):
        print("fidelity is {0}".format(self.fidelity))
        plt.figure(figsize=(5,5))
        xfit = np.linspace(0,self.standard_rb.length_list[-1],1001)
        srb_yfit = exp_decay(xfit,self.standard_rb.a,self.standard_rb.b,self.standard_rb.p)
        irb_yfit = exp_decay(xfit,self.interleaved_rb.a,self.interleaved_rb.b,self.interleaved_rb.p)
        plt.plot(xfit,srb_yfit,'r-', label="Standard")
        plt.plot(xfit,irb_yfit,'b-', label="Interleaved")
        plt.errorbar(x=self.standard_rb.length_list, y=self.standard_rb.pauli_ave, yerr=self.standard_rb.pauli_std, fmt='k.')
        plt.errorbar(x=self.interleaved_rb.length_list, y=self.interleaved_rb.pauli_ave, yerr=self.interleaved_rb.pauli_std, fmt='k.')
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.legend()
        plt.show()