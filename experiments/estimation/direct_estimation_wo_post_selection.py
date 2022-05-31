import copy
from ...objects import Job, JobTable, Report

class DirectEstimation:
    def __init__(
        self,
        ansatz,
        circuit,
        qubit_index,
        spam_condition_list,
        ):

        self.name = "DirectEstimation"
        
        self.job_table  = JobTable(name=self.name)
        for condition in spam_condition_list:
            
            cir = copy.deepcopy(circuit)
            for i, (pauli, index) in enumerate(zip(condition["prep_pauli"], condition["prep_index"])):
                cir.prep_init(pauli, index, qubit_index[i])
            ansatz(cir)
            for i, pauli in enumerate(condition["meas_pauli"]):
                cir.meas_axis(pauli, qubit_index[i])
            cir.qtrigger(qubit_index)
            cir.measurement_all()
            
            condition["sequence"] = cir
            self.job_table.submit(Job(condition))

    def execute(self, take_data):
        take_data(self.job_table)

    def make_data_table(self):
        self.data_table = {}
        for job in self.job_table.table:
            if (job.prep_pauli, job.meas_pauli) not in self.data_table.keys():
                self.data_table[(job.prep_pauli, job.meas_pauli)] = {}
            self.data_table[(job.prep_pauli, job.meas_pauli)][job.prep_index] = job.result

    def reset(self):
        self.job_table.reset()
        self.data_table = {}
        self.report = None

class MitigatedDirectEstimation(DirectEstimation):
    def __init__(
        self,
        ansatz,
        circuits,
        spam_condition_list,
        ):

        self.name = "MitigatedDirectEstimation"
        self.job_tables = {}
        for key, circuit in circuits.items():
            job_table  = JobTable()
            for spam_condition in spam_condition_list:
                circuit._reset()
                for i, (pauli, index) in enumerate(zip(spam_condition["prep_pauli"], spam_condition["prep_index"])):
                    circuit.state_preparation(pauli, index, circuit.qubit_name_list[i])
#                 cir.call(ansatz)
                ansatz(cir)
                for i, pauli in enumerate(spam_condition["meas_pauli"]):
                    circuit.measurement(pauli, circuit.qubit_name_list[i])
                spam_condition["sequence"] = copy.deepcopy(circuit.base.sequence)
                job_table.submit(Job(spam_condition))
            self.job_tables[key] = job_table

    def execute(self, take_data):
        for job_table in self.job_tables.values():
            take_data(job_table)

    def make_data_table(self):
        self.data_tables = {}
        for key, job_table in self.job_tables.items():
            data_table = {}
            for job in job_table.table:
                if (job.prep_pauli, job.meas_pauli) not in data_table.keys():
                    data_table[(job.prep_pauli, job.meas_pauli)] = {}
                data_table[(job.prep_pauli, job.meas_pauli)][job.prep_index] = job.result
            self.data_tables[key] = data_table

    def reset(self):
        for job_table in self.job_tables.values():
            job_table.reset()
        for data_table in self.data_tables.values():
            data_table = {}
        self.report = None
