import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report

class QuantumProcessTomography:
    def __init__(
        self,
        number_of_qubit
        ):

        self.name = "QuantumProcessTomography"
        self.number_of_qubit = number_of_qubit

    def set_circuit(self, circuits, qubit_index):
        self.circuits = circuits
        self.qubit_index = qubit_index

    def prepare(self, ansatz):
        spam_condition_list = []
        
        spam_1q = list(itertools.product(["X","Y","Z"], ["0","1"], ["X", "Y", "Z"]))
        spam_nq = list(itertools.product(spam_1q, repeat=self.number_of_qubit))
        for spam in spam_nq:
            prep_pauli = ""
            prep_index = ""
            meas_pauli = ""
            for i in spam:
                prep_pauli += i[0]
                prep_index += i[1]
                meas_pauli += i[2]
            spam_condition_list.append(
                {
                    "prep_pauli" : prep_pauli,
                    "meas_pauli" : meas_pauli,
                    "prep_index" : prep_index,
                }
            )
                
        self.des = {}
        for key, circuit in self.circuits.items():
            self.des[key] = DirectEstimation(ansatz, circuit, self.qubit_index, spam_condition_list)

    def execute(self, take_data):
        for de in self.des.values():
            de.execute(take_data)

    def analyze(self):
        for de in self.des.values():
            de.make_data_table()

        self.report = Report(name="quantum_process_tomography")
        for key, val in self.des.items():
            self.report.add_information(f"data table {key}", val.data_table)

    def visualize(self):
        print(self.report.name)