import numpy as np
from rotated_bridge.qubit import Qubit
from rotated_bridge.gates import Gate, Reset, Hadamard, Measurement, Cnot


class Hexagon():
    def __init__(self, position: np.ndarray, qubits: np.ndarray , index: int) -> None:
        self.position = position
        self.qubits = qubits
        self.index = index

    def __str__(self) -> str:
        """Returns string representing hexagon"""
        string = ""
        for row in self.qubits:
            for q in row:
                if isinstance(q, Qubit):
                    index = q.index
                    string += f"{index:03},"
                else:
                    string += "|-|,"
            string += "\n"

        return string
    
    def __repr__(self) -> str:
        """Returns string representing hexagon"""
        return f"Hexagon(pos{self.position}, index={self.index})"
    
    def fold_top(self, start_time: int, cycle=0):
        qubits = {
            "a" : self.qubits[(2, 2-2*cycle)],
            "b" : self.qubits[(1, 2-2*cycle)],
            "c" : self.qubits[(0, 2-2*cycle)],
            "d" : self.qubits[(0, 1)],
            "e" : self.qubits[(0, 0+2*cycle)],
            "f" : self.qubits[(1, 0+2*cycle)],
            "g" : self.qubits[(2, 0+2*cycle)]
        }

        schedule = []

        clock = 0
        r = Reset(target=qubits["d"], time=clock + start_time)
        if(r.isvalid()): schedule.append(r)
        
        clock += 1 # 1
        c1 = Cnot(control=qubits["b"], target=qubits["a"], time=clock + start_time)
        c2 = Cnot(control=qubits["e"], target=qubits["f"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)

        clock += 1 # 2
        c1 = Cnot(control=qubits["a"], target=qubits["b"], time=clock + start_time)
        c2 = Cnot(control=qubits["f"], target=qubits["e"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)

        clock += 1 # 3
        c1 = Cnot(control=qubits["b"], target=qubits["c"], time=clock + start_time)
        c2 = Cnot(control=qubits["g"], target=qubits["f"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)
        
        clock += 1 # 4
        c1 = Cnot(control=qubits["a"], target=qubits["b"], time=clock + start_time)
        c2 = Cnot(control=qubits["c"], target=qubits["d"], time=clock + start_time)
        c3 = Cnot(control=qubits["f"], target=qubits["e"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)
        if(c3.isvalid()): schedule.append(c3)

        clock += 1 # 5
        c1 = Cnot(control=qubits["e"], target=qubits["d"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
 
        clock += 1 # 6
        c1 = Cnot(control=qubits["a"], target=qubits["b"], time=clock + start_time)
        c2 = Cnot(control=qubits["f"], target=qubits["e"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)

        clock += 1 # 7
        c1 = Cnot(control=qubits["b"], target=qubits["c"], time=clock + start_time)
        c2 = Cnot(control=qubits["g"], target=qubits["f"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)

        clock += 1 # 8
        c1 = Cnot(control=qubits["a"], target=qubits["b"], time=clock + start_time)
        c2 = Cnot(control=qubits["f"], target=qubits["e"], time=clock + start_time)
        if c1.isvalid(): schedule.append(c1)
        if c2.isvalid(): schedule.append(c2)

        clock += 1 # 9
        c1 = Cnot(control=qubits["b"], target=qubits["a"], time=clock + start_time)
        c2 = Cnot(control=qubits["e"], target=qubits["f"], time=clock + start_time)
        if c1.isvalid(): 
            schedule.append(c1)
        if c2.isvalid(): 
            schedule.append(c2)

        clock += 1 # 10
        m = Measurement(target=qubits["d"], time=clock + start_time, xz="z") 
        if(m.isvalid()): schedule.append(m)

        return schedule


    def fold_bottom(self, start_time: int, cycle=0) -> list[Gate]:
        """
        initialises list of appropriate gates to simulate folding x stabilisers into two body
        """
        qubits = {
            "e" : self.qubits[(2, 2-2*cycle)],
            "f" : self.qubits[(3, 2-2*cycle)],
            "g" : self.qubits[(4, 2-2*cycle)],
            "h" : self.qubits[(4, 1)],
            "i" : self.qubits[(4, 0+2*cycle)],
            "j" : self.qubits[(3, 0+2*cycle)],
            "k" : self.qubits[(2, 0+2*cycle)]
        }

        schedule = []

        clock = 0
        r = Reset(target=qubits["h"], time=clock+start_time)
        if(r.isvalid()): schedule.append(r)

        clock += 1 # 1
        c1 = Cnot(control= qubits["e"], target=qubits["f"], time=clock + start_time)
        c2 = Cnot(control= qubits["j"], target=qubits["i"], time=clock + start_time)
        h = Hadamard(target=qubits["h"], time=clock+start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)
        if(h.isvalid()): schedule.append(h)

        clock += 1 # 2
        c1 = Cnot(control= qubits["f"], target=qubits["e"], time=clock + start_time)
        c2 = Cnot(control= qubits["i"], target=qubits["j"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)

        clock += 1 # 3
        c1 = Cnot(control= qubits["g"], target=qubits["f"], time=clock + start_time)
        c2 = Cnot(control= qubits["j"], target=qubits["k"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)

        clock += 1 # 4
        c1 = Cnot(control= qubits["f"], target=qubits["e"], time=clock + start_time)
        c2 = Cnot(control= qubits["h"], target=qubits["g"], time=clock + start_time)
        c3 = Cnot(control= qubits["i"], target=qubits["j"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)
        if(c3.isvalid()): schedule.append(c3)
        
        clock += 1 # 5
        c1 = Cnot(control= qubits["h"], target=qubits["i"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        
        clock += 1 # 6
        c1 = Cnot(control= qubits["f"], target=qubits["e"], time=clock + start_time)
        c2 = Cnot(control= qubits["i"], target=qubits["j"], time=clock + start_time)
        h = Hadamard(target=qubits["h"], time=clock+start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)
        if(h.isvalid()): schedule.append(h)

        clock += 1 # 7
        c1 = Cnot(control= qubits["g"], target=qubits["f"], time=clock + start_time)
        c2 = Cnot(control= qubits["j"], target=qubits["k"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)

        clock += 1 # 8
        c1 = Cnot(control= qubits["f"], target=qubits["e"], time=clock + start_time)
        c2 = Cnot(control= qubits["i"], target=qubits["j"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)

        clock += 1 # 9
        c1 = Cnot(control= qubits["e"], target=qubits["f"], time=clock + start_time)
        c2 = Cnot(control= qubits["j"], target=qubits["i"], time=clock + start_time)
        if(c1.isvalid()): schedule.append(c1)
        if(c2.isvalid()): schedule.append(c2)

        clock += 1 # 10
        m = Measurement(target=qubits["h"], time=clock + start_time, xz="x") 
        if(m.isvalid()): schedule.append(m)

        return schedule
    
    def get_final_detectors(self, measurements, m_count) -> str:
        top = {
            1 : self.qubits[0,0], 
            2 : self.qubits[0,2], 
            3 : self.qubits[2,2], 
            4 : self.qubits[2,0],
            -1: self.qubits[0,1]
            }

        string = "DETECTOR"
        for q in top.values():
            x,y,_=q
            m = measurements[y,x][-1]
            if(m.xz == 'x'):
                m = measurements[y,x][-2]
            _,_,_,count = m
            i = m_count-count
            string += f" rec[{-i}]"
        
        return string + "\n"


        


        