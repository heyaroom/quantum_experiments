import itertools
import numpy as np
from .direct_estimation import DirectEstimation
from ...objects import Report

class QuantumStateTomography:
    def __init__(
        self,
        number_of_qubit
        ):

        self.name = "QuantumStateTomography"
        self.number_of_qubit = number_of_qubit

    def set_circuit(self, circuits, qubit_index):
        self.circuits = circuits
        self.qubit_index = qubit_index

    def prepare(self, ansatz):
        spam_condition_list = []
        
        spam_1q = ["X", "Y", "Z"]
        spam_nq = list(itertools.product(spam_1q, repeat=self.number_of_qubit))
        for spam in spam_nq:
            meas_pauli = ""
            for i in spam:
                meas_pauli += i[0]
            spam_condition_list.append(
                {
                    "prep_pauli" : "I"*self.number_of_qubit,
                    "meas_pauli" : meas_pauli,
                    "prep_index" : "0"*self.number_of_qubit,
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

        self.report = Report(name="quantum_state_tomography")
        for key, val in self.des.items():
            self.report.add_information(f"data table {key}", val.data_table)

    def visualize(self):
        print(self.report.name)