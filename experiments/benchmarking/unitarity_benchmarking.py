import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ...objects import Report, Job, JobTable

def exp_decay(x,a,b,p):
    y = a*p**x + b
    return y

def double_exp_decay(x,p1,p2):
    y = 1/3.*p1**x + 2/3.*p2**x
    return y

class UnitarityBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_index,
        group,
        sequence_list,
        seed = 0,
        interleaved = None,
        ):

        self.name               = "UnitarityBenchmarking"
        self.number_of_qubit    = group.num_qubit
        self.qubit_index        = qubit_index
        self.seed               = seed
        self.interleaved        = interleaved
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        self.random_index       = np.array(sequence_list).T[1].tolist()[0]
        self.report             = Report(name="unitarity_randomized_benchmarking")
        
        self.job_table = JobTable(name=self.name)
        for (length, random, shot) in self.sequence_list:
            
            ## generate gate_array ##
            sequence_array = []
            rand_gate_array = group.sample(random*(length-1), seed=self.seed)
            rand_gate_array = rand_gate_array.reshape(random, length-1, 2**self.number_of_qubit, 2**self.number_of_qubit)
            for rand_gates in rand_gate_array:
                gate_array = []
                gate = np.identity(2**self.number_of_qubit)
                for rand in rand_gates:
                    gate_array.append(rand)
                    gate = rand@gate
                    if self.interleaved is not None:
                        gate = self.interleaved["gate"]@gate
                sequence_array.append(gate_array)
                    
            ## apply experiment ##
            for gate_array in sequence_array:
                cir = copy.deepcopy(circuit)
                for pos, gate in enumerate(gate_array):
                    if int(np.log2(gate.shape[0])) == 1:
                        cir.su2(gate, target=self.qubit_index[0])
                    if int(np.log2(gate.shape[0])) == 2:
                        cir.su4(gate, control=self.qubit_index[0], target=self.qubit_index[1])
                    if interleaved is not None:
                        if pos != len(gate_array) - 1:
                            cir.qtrigger(self.qubit_index)
                            cir.call(self.interleaved["ansatz"])
                            cir.qtrigger(self.qubit_index)
                cir.qtrigger(self.qubit_index)
                cir.measurement_all()

                ## job submition ##
                condition = {
                    "length"         : length,
                    "gate_array"     : gate_array,
                    "shot"           : shot,
                    "sequence"       : cir,
                }
                self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)
        
    def tmp_analyze(self):
        self.hist_table = {}
        for length in self.length_list:
            self.hist_table[length] = []
        for job in self.job_table.table:
            self.hist_table[job.length].append(job.result)
            
        self.report.add_information("hist table", self.hist_table)
        self.report.add_information("sequence", self.sequence_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit index", self.qubit_index)

    def analyze(self):
        
        if self.number_of_qubit == 1:
            pauli = np.array([job.result["0"] - job.result["1"] for job in self.job_table.table])
        elif self.number_of_qubit == 2:
            pauli = np.array([job.result["00"] + job.result["11"] - job.result["01"] - job.result["10"] for job in self.job_table.table])
            
        pauli = pauli.reshape(len(self.length_list), self.random_index)
        self.pauli = pauli

        self.b = np.var(pauli, axis=1) #*(4**self.number_of_qubit - 1)

        popt, pcov = curve_fit(exp_decay,self.length_list,self.b,p0=[self.b[0],0,0.99])
        self.b_fit_param = popt
        self.unitarity = popt[2]
        
        self.report.add_information("fit params for unitarity", self.b_fit_param)
        self.report.add_information("unitarity", self.unitarity)
        self.report.add_information("pauli", self.pauli)

    def visualize(self):
        popt = self.b_fit_param
        bfit = exp_decay(self.length_list, popt[0], popt[1], popt[2])
        plt.figure(figsize=(5,5))
        plt.plot(self.length_list, self.b, "b.", label=f"unitarity {self.unitarity}")
        plt.plot(self.length_list, bfit, "b-")
        plt.axhline(0, color="black", linestyle="--")
        plt.ylim(-1,1)
        plt.legend()
        plt.show()