import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ...objects import Report, Job, JobTable

def exp_decay(x,a,b,p):
    y = a*p**x + b
    return y

def pre_process(sequences):
    I = np.identity(2, dtype=np.complex128)
    
    random = sequences.shape[0]
    length = sequences.shape[1] + 1
    number_of_qubit = sequences.shape[2]

    post_sequences = np.zeros([random, length, number_of_qubit, 2, 2], dtype=np.complex128)
    post_sequences[:, :-1] = sequences

    for sequence, post_sequence in zip(sequences, post_sequences):
        gate_prods = [I]*number_of_qubit
        for gate_list in sequence:
            for i in range(number_of_qubit):
                gate_prods[i] = gate_list[i]@gate_prods[i]
                
        global_inverse = [gate.T.conj() for gate in gate_prods]
        post_sequence[-1] = global_inverse
        
    return post_sequences

class SimultaneousRandomizedBenchmarking:
    def __init__(
        self,
        circuit,
        qubit_list,
        group,
        sequence_list,
        seed = 0,
        ):

        self.name               = "SimultaneousRandomizedBenchmarking"
        self.number_of_qubit    = 1
        self.qubit_list         = qubit_list
        self.seed               = seed
        self.sequence_list      = sequence_list
        self.length_list        = np.array(sequence_list).T[0].tolist()
        self.random_list      = np.array(sequence_list).T[1].tolist()
        self.report             = Report(name=self.name)
        self.job_table          = JobTable(name=self.name)

        for (length, random, shot) in self.sequence_list:
            if length == 0:
                for idx in range(random):
                    cir = copy.deepcopy(circuit)
                    cir.qtrigger(self.qubit_list)
                    cir.measurement_all()

                    condition = {
                        "length"      : length,
                        "shot"        : shot,
                        "sequence"    : cir,
                    }
                    self.job_table.submit(Job(condition))

            else:
                rand_sequences = group.sample(random*(length-1)*len(self.qubit_list), seed=self.seed)
                rand_sequences = rand_sequences.reshape(random, length-1, len(self.qubit_list), 2, 2)
                rand_sequences = pre_process(rand_sequences)

                for tmp_sequence in rand_sequences:
                    cir = copy.deepcopy(circuit)
                    for tmp_gates in tmp_sequence:
                        for qubit_index, tmp_gate in zip(self.qubit_list, tmp_gates):
                            cir.su2(tmp_gate, target=qubit_index)
                        cir.qtrigger(self.qubit_list)
                    cir.measurement_all()

                    condition = {
                        "length"      : length,
                        "shot"        : shot,
                        "sequence"    : cir,
                    }
                    self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)
        
    def analyze(self):

        result = [job.result for job in self.job_table.table]
        
        population_list = []
        for res in result:
            population = [0]*len(self.qubit_list)
            for key, val in res.items():
                for i in range(len(self.qubit_list)):
                    if key[i] == "0":
                        population[i] += val
            population_list.append(population)
        population_list = np.array(population_list).reshape(len(self.length_list), self.random_list[0], len(self.qubit_list))

        self.ave_list = np.mean(population_list, axis=1).T
        self.std_list = np.std(population_list, axis=1).T
        
        self.fit_params = []
        for ave in self.ave_list:
            a0 = 0.5
            b0 = 0.5

            _x = np.array(self.length_list)
            _f = (ave - b0)/a0
            _x = _x[np.where(_f>0)]
            _f = _f[np.where(_f>0)]
            _y = np.log(_f)
            p0 = np.exp(np.mean(np.gradient(_y,_x)))
            popt, pcov = curve_fit(exp_decay, self.length_list, ave, p0=[a0,b0,p0])
            self.fit_params.append(popt)
        
        self.fidelities = [0.5*(1+p) for (a,b,p) in self.fit_params]

        self.report.add_information("average gate fidelties", self.fidelities)
        self.report.add_information("fit params : a, b, p", self.fit_params)
        self.report.add_information("population : average", self.ave_list)
        self.report.add_information("population : standard deviation", self.std_list)
        self.report.add_information("seed", self.seed)
        self.report.add_information("qubit list", self.qubit_list)

    def visualize(self):
        for qubit_index, fidelity in zip(self.qubit_list, self.fidelities):
            print(f"qubit : {qubit_index}, Clifford fidelity : {fidelity}")
            print(f"qubit : {qubit_index}, RX90 fidelity : {0.5*(1+(2*fidelity-1)**0.5)}")
            
        plt.figure(figsize=(5,5))
        xfit = np.linspace(0,self.length_list[-1],1001)
        for idx, (a,b,p) in enumerate(self.fit_params):
            yfit = exp_decay(xfit, a, b, p)
            plt.plot(xfit, yfit, color=f"C{idx}")
        for idx, (ave, std) in enumerate(zip(self.ave_list, self.std_list)):
            plt.errorbar(x=self.length_list, y=ave, yerr=std, fmt=".", color=f"C{idx}", label=f"Q{self.qubit_list[idx]}")
        plt.axhline(0.5, color="black", linestyle="--")
        plt.xlabel('Sequence length')
        plt.ylabel("Population")
        plt.ylim(-0.1,1.1)
        plt.show()

    def reset(self):
        self.job_table.reset()
        self.data_table = {}
        self.report = None