from rotated_bridge.gates import Gate, Cnot, Hadamard, Measurement, Reset
import numpy as np


def clean_gate_list(ls : list[Gate]) -> list[Gate]:
    clean_ls = []
    for gate in ls:
        if gate not in clean_ls:
            clean_ls.append(gate)

    min_time = int(1e100)
    max_time = -int(1e100)

    for gate in ls:
        if gate.time > max_time:
             max_time = int(gate.time)
        if gate.time < min_time:
             min_time = int(gate.time)

    ordered = []
    for t in range(min_time, max_time+1):
        for gate in clean_ls:
            if(gate.time == t):
                ordered.append(gate)
    
    return ordered
    

def gate_list_to_string(ls : list[Gate]) -> str:
    min_time = int(1e100)
    max_time = -int(1e100)

    for gate in ls:
        if gate.time > max_time:
             max_time = int(gate.time)
        if gate.time < min_time:
             min_time = int(gate.time)

    string = ""
    
    for t in range(min_time, max_time+1):

        reset = "R"
        rx_error = "X_ERROR(0.02)"
        cnot = "CNOT"
        c_depolarise = "DEPOLARIZE1(0.02)"
        measure = "M"
        mx_error = "X_ERROR(0.02)"
        hadamard = "H"
        h_depolarise = "DEPOLARIZE1(0.02)"

        for gate in ls:

            if(gate.time == t):
               if isinstance(gate, Reset):
                    reset += f" {gate.get_target_index()}"
                    rx_error += f" {gate.get_target_index()}"
               elif isinstance(gate, Cnot):
                    cnot += f" {gate.get_control_index()} {gate.get_target_index()}"
                    c_depolarise += f" {gate.get_control_index()} {gate.get_target_index()}"
               elif isinstance(gate, Measurement):
                    measure += f" {gate.get_target_index()}"
                    mx_error += f" {gate.get_target_index()}"
               elif isinstance(gate, Hadamard):
                    hadamard += f" {gate.get_target_index()}"
                    h_depolarise += f" {gate.get_target_index()}"
               else:
                   raise ValueError
                

        if(len(hadamard) > 1):
             string += hadamard + "\n"
             string += h_depolarise+ "\n"
        if(len(reset) > 1):
             string += reset + "\n"
             string += rx_error+ "\n"
        if(len(cnot) > 4):
             string += cnot + "\n"
             string += c_depolarise+ "\n"
        if(len(measure) > 1):
             string += mx_error+ "\n"
             string += measure + "\n"
        string +=  "TICK\n"

    return string


def find_measurement(target: Measurement, ls: list[Gate]) -> list[int]:

    m_ls = [m for m in ls if isinstance(m, Measurement)]

    indices = []
    for i, m in enumerate(m_ls):
        if(m == target):
            indices.append(-(len(m_ls) - i))

    if(len(indices) > 0):
         return indices
    else:
         raise ValueError