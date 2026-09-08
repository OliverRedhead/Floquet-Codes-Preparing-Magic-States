from rotated_bridge.qubit import Qubit, Data, Bridge, Syndrome, SyData
import numpy as np
from rotated_bridge.gates import Hadamard, Cnot, Measurement, Reset, Gate
from rotated_bridge.helpers import clean_gate_list, find_measurement
from rotated_bridge.hexagon import Hexagon
import stim


class Surface():

    def __init__(self, dist, reps=1, errors=[0.0002,0.0041,.042,.012,.075]) -> None:
        
        if dist%2 == 0:
            raise NotImplementedError
        
        self.clock = 0
        self.reps = reps
        self.dist = dist
        self.measurements = []

        self.errors = {
            "single_qubit" : errors[0],
            "two_qubit" : errors[1],
            "measurement" : errors[2],
            "idle" : errors[3],
            "reset" : errors[4]
        }

        self.initialise_qubits()
        self.set_hexagons()
        self.initialise_measurements()

    def set_hexagons(self):
        """
        set hexagons on self.qubits
        """
        n = self.dist*2
        canvas = np.empty(shape=(n+6, n+3), dtype=Qubit)
        canvas[3:-2, 1:-1] = self.qubits
        hexagons = []

        for x in range(2, n-1, 2):
            ystart = 4
            if(x%4==0):
                ystart = 2

            for y in range(ystart, n+1, 4):
                
                if (x,y) == (0,n): continue

                x_hex = x
                y_hex = y+1
                position = np.array([x,y])

                hex_qubits = canvas[y_hex:y_hex+5,x_hex:x_hex+3]
                h = Hexagon(position=position, qubits=hex_qubits, index=len(hexagons))
                hexagons.append(h)            
        self.hexagons = hexagons 
        

    def initialise_qubits(self) -> None:
        """
        initialise dist+1 x dist+1 grid of qubits in heavy hex lattice
        """
        n = 2*self.dist
        grid = np.empty(shape=(n+1, n+1), dtype=Qubit)
        
        count = 0
        for x in range(0, n+1):
            # this is for left and right edges (syndrome cases)
            ystart = 0 
            if x == 0 or x == n: 
                ystart = 3

            for y in range(ystart, n+1):

                if (x,y) == (1,0) or (x,y) == (1,1):
                    continue

                position = np.array([x,y])

                if(x%2 == 0): 
                    if((x+y)%4 == 0):
                        grid[y,x] = Syndrome(position=position, index=count)
                        count += 1
                    continue

                elif(y%2 == 0): # Data
                    if(y == 0):
                        grid[y,x] = SyData(position=position, index=count)
                    else: 
                        grid[y,x] = Data(position=position, index=count)
                else: # bridge
                    grid[y,x] = Bridge(position=position, index=count)
                
                count += 1

        self.qubits = grid

    def initialise_measurements(self) -> None:
        """
        returns intialised nupmy array with list[measurement] data types. 
        This will hold measurements so I can access their coordinates more easily
        """
        shape = self.qubits.shape
        arr = np.empty(shape, dtype=object)

        for i in range(shape[0]):
            for j in range(shape[1]):
                arr[i, j] = []

        self.measurements = arr
        self.m_count = 0

    def set_measurements(self, gates_ls: list[Gate]) -> None:
        """
        sets all measurement instances iin gates_ls onto class attribute measurements
        """
        for gate in gates_ls:
            if isinstance(gate, Measurement):
                gate.set_count(self.m_count)
                self.m_count+=1
                x, y = gate.target.position[0], gate.target.position[1]
                self.measurements[y,x].append(gate)             # type: ignore

    def get_circuit(self):
        string = ""
        for row in self.qubits:
            for q in row:
                if isinstance(q, Qubit):
                    string += str(q) + "\n"
        
        gates = self.reset_data()
        string += self.measure_stabilisers1(gates)
        string += self.measure_stabilisers2()
        string += self.intitial_detectors()

        string += "REPEAT " + str(self.reps) + "{\n"

        string += self.measure_stabilisers1()
        string += self.measure_stabilisers2()
        string += self.get_detectors()

        string += "\n}\n"

        string += self.measure_data()
        string += self.final_detectors()
        string += self.observable()


        return string
    
    def reset_data(self) -> list:
        gates = []
        for row in self.qubits:
            for q in row:
                if isinstance(q, Data) and not isinstance(q,SyData):
                    gates.append(Reset(target=q, time=0))
        return gates

    def measure_stabilisers1(self, gates_ls=None) -> str:
        if gates_ls is None:
            gates_ls = []
        for h in self.hexagons:
            if(h.position[0]%4 == 0):
                gates_ls += h.fold_top(start_time=self.clock,cycle=0)
            else:
                gates_ls += h.fold_bottom(start_time=self.clock,cycle=0)

        for q in self.qubits[:,-1]: # right edge
            if isinstance(q, Qubit):
                x,y,i=q
                gates_ls.append(Reset(target=q,time=self.clock))
                gates_ls.append(Hadamard(target=q,time=self.clock+1))
                gates_ls.append(Cnot(control=q,target=self.qubits[y,x-1],time=self.clock+5))
                gates_ls.append(Hadamard(target=q,time=self.clock+6))
                gates_ls.append(Measurement(target=q,time=self.clock+10,xz="x"))

        for q in self.qubits[0,:]:
            if isinstance(q,SyData):
                x,y,i = q
                gates_ls.append(Reset(target=q,time=self.clock))

        gates_ls = clean_gate_list(gates_ls)
        self.set_measurements(gates_ls)

        start_time = self.clock
        end_time = self.clock + 11
        string = ""
        for t in range(start_time, end_time):
            
            hadamard = ""
            reset = ""
            measurement = ""
            cnot = ""

            for gate in gates_ls:
                if gate.time == t:
                    if isinstance(gate, Measurement):
                        measurement += f" {gate.get_target_index()}"
                    elif isinstance(gate, Hadamard):
                        hadamard += f" {gate.get_target_index()}"
                    elif isinstance(gate, Reset):
                        reset += f" {gate.get_target_index()}"
                    elif isinstance(gate, Cnot):
                        cnot += f" {gate.get_control_index()} {gate.get_target_index()}"

            if(len(measurement) > 0):
                string += "M" + measurement + "\n"
                string += f"X_ERROR({self.errors['measurement']})" + measurement + "\n"
                string += f"DEPOLARIZE1({self.errors['idle']})" + self.all_other_qubits(measurement) + "\n"
                string += "\nTICK\n"
                continue
            if len(reset) > 1:
                string += "R" + reset + "\n"
                string += f"X_ERROR({self.errors['reset']})" + reset + "\n"
                string += f"DEPOLARIZE1({self.errors['idle']})" + self.all_other_qubits(reset) + "\n"
                string += "\nTICK\n"
                continue
            
            if len(hadamard) > 1:
                string += "H" + hadamard + "\n"
                string += f"DEPOLARIZE1({self.errors['single_qubit']})" + hadamard + "\n"
            if len(cnot) > 1:
                string += "CNOT" + cnot + "\n"
                string += f"DEPOLARIZE2({self.errors['two_qubit']})" + cnot + "\n"
            
            string += "\nTICK\n"
    
        return string
    
    def measure_stabilisers2(self) -> str:
        gates_ls = []

        for h in self.hexagons:
            if(h.position[0]%4 == 2):
                gates_ls += h.fold_top(start_time=self.clock,cycle=1)
            else:
                gates_ls += h.fold_bottom(start_time=self.clock,cycle=1)

        for q in self.qubits[:,0]: # left edge
            if isinstance(q, Qubit):
                x,y,i=q
                gates_ls.append(Reset(target=q,time=self.clock))
                gates_ls.append(Hadamard(target=q,time=self.clock+1))
                gates_ls.append(Cnot(control=q,target=self.qubits[y,x+1],time=self.clock+5))
                gates_ls.append(Hadamard(target=q,time=self.clock+6))
                gates_ls.append(Measurement(target=q,time=self.clock+10,xz="x"))

        # top and bottom edge cases
        for q in self.qubits[0,:]:
            if isinstance(q,Data):
                gates_ls.append(Measurement(target=q,time=self.clock+10,xz="z"))

        for q in self.qubits[-1,:-1]:
            if isinstance(q,Syndrome):
                x,y,i=q
                gates_ls.append(Reset(target=q,time=self.clock))
                gates_ls.append(Cnot(control=self.qubits[y,x-1],target=q,time=self.clock+5))
                gates_ls.append(Cnot(control=self.qubits[y,x+1],target=q,time=self.clock+6))
                gates_ls.append(Measurement(target=q,time=self.clock+10,xz="z"))

        gates_ls = clean_gate_list(gates_ls)
        self.set_measurements(gates_ls)

        start_time = self.clock
        end_time = self.clock + 11
        string = ""
        for t in range(start_time, end_time):
            
            hadamard = ""
            reset = ""
            measurement = ""
            cnot = ""

            for gate in gates_ls:
                if gate.time == t:
                    if isinstance(gate, Measurement):
                        measurement += f" {gate.get_target_index()}"
                    elif isinstance(gate, Hadamard):
                        hadamard += f" {gate.get_target_index()}"
                    elif isinstance(gate, Reset):
                        reset += f" {gate.get_target_index()}"
                    elif isinstance(gate, Cnot):
                        cnot += f" {gate.get_control_index()} {gate.get_target_index()}"

            if(len(measurement) > 0):
                string += "M" + measurement + "\n"
                string += f"X_ERROR({self.errors['measurement']})" + measurement + "\n"
                string += f"DEPOLARIZE1({self.errors['idle']})" + self.all_other_qubits(measurement) + "\n"
                string += "\nTICK\n"
                continue
            if len(reset) > 1:
                string += "R" + reset + "\n"
                string += f"X_ERROR({self.errors['reset']})" + reset + "\n"
                string += f"DEPOLARIZE1({self.errors['idle']})" + self.all_other_qubits(reset) + "\n"
                string += "\nTICK\n"
                continue
    
            if len(hadamard) > 1:
                string += "H" + hadamard + "\n"
                string += f"DEPOLARIZE1({self.errors['single_qubit']})" + hadamard + "\n"
            if len(cnot) > 1:
                string += "CNOT" + cnot + "\n"
                string += f"DEPOLARIZE2({self.errors['two_qubit']})" + cnot + "\n"
            
            string += "\nTICK\n"
    
        return string
    
    def all_other_qubits(self, measurement_string) -> str:
        string = ""
        indices_ls = [int(x) for x in measurement_string.split()]
        for row in self.qubits:
            for q in row:
                if isinstance(q, Qubit) and q.index not in indices_ls:
                    string += f" {q.index}"
        return string

    def intitial_detectors(self) -> str:
        string = ""
        for row in self.measurements[1:]: # skip top edge - will have to do this differently
            for m_ls in row:
                if len(m_ls) == 1 and m_ls[0].xz == 'z':
                    m = m_ls[0]    
                elif len(m_ls) == 2 and m_ls[0].xz == 'z':
                    m = m_ls[0]
                elif len(m_ls) == 2 and m_ls[1].xz == 'z':
                    m = m_ls[1]
                else:
                    continue

                i = self.m_count - m.count 
                x,y,t = m.position[0], m.position[1], m.time
                string += f"DETECTOR({x}, {y}, {t}) rec[{-i}]\n"

        ms = []
        for m_ls in self.measurements[0]:
            if len(m_ls) > 0:
                ms.append(m_ls[-1])
        
        for i in range(0, len(ms), 3):
            m1 = ms[i+1]
            m2 = ms[i]
            m3 = ms[i+2]

            j = self.m_count - m1.count 
            x,y,t = m1.position[0], m1.position[1], m1.time
            string += f"DETECTOR({x}, {y}, {t}) rec[{-j}]"
            j = self.m_count - m2.count 
            string += f" rec[{-j}]"
            j = self.m_count - m3.count 
            string += f" rec[{-j}]\n"
            
        return string

    def get_detectors(self) -> str:
        string = ""
        measurements = np.array(self.measurements)
        for row in self.measurements[1:]:
            for m_ls in row:
                if len(m_ls) == 0:
                    continue
                
                if len(m_ls) == 1:
                    m = m_ls[0]
                    x,y,t,count = m
                    i = self.m_count - count # TODO check this
                    string += f"DETECTOR({x},{y},{t}) rec[-{i}]\n"
                    continue

                m1 = m_ls[-1]
                m2 = m_ls[-2]
                if(m1.xz == m2.xz):
                    x1,y1,t1,count1 = m1
                    x2,y2,t2,count2 = m2
                    i1 = self.m_count - count1 # TODO check this
                    i2 = self.m_count - count2 # TODO check this
                    string += f"DETECTOR({x1},{y1},{t1}) rec[-{i1}] rec[-{i2}]\n"
                else:
                    m1a = m1
                    m1b = m_ls[-3]

                    m2a = m2
                    m2b = m_ls[-4]

                    x1,y1,t1,count1 = m1a
                    x2,y2,t2,count2 = m1b
                    i1 = self.m_count - count1 # TODO check this
                    i2 = self.m_count - count2 # TODO check this

                    string += f"DETECTOR({x1},{y1},{t1}) rec[-{i1}] rec[-{i2}]\n"

                    x1,y1,t1,count1 = m2a
                    x2,y2,t2,count2 = m2b
                    i1 = self.m_count - count1 # TODO check this
                    i2 = self.m_count - count2 # TODO check this

                    string += f"DETECTOR({x1},{y1},{t1}) rec[-{i1}] rec[-{i2}]\n"

        for m_ls in self.measurements[0]: # 6 body at top
            if len(m_ls) == 0 or isinstance(m_ls[0].target, SyData):
                continue

            m_syndrome1 = m_ls[-1]
            m_syndrome2 = m_ls[-2]
            x,y,t,_ = m_syndrome1
            m_sydata1_1 = measurements[0,x-1][-1]
            m_sydata1_2 = measurements[0,x-1][-2]
            m_sydata2_1 = measurements[0,x+1][-1]
            m_sydata2_2 = measurements[0,x+1][-2]

            m_list = [m_syndrome1,m_syndrome2, m_sydata1_1,m_sydata1_2,m_sydata2_1,m_sydata2_2]

            indices = [ self.m_count-m.count for m in m_list ]
            string += f"DETECTOR({x},{y},{t})"
            for i in indices:
                string += f" rec[{-i}]"
            string += "\n"

        return string    
    
    # TODO measure 6 body of N syndrome N between rounds

    def measure_data(self):
        m_ls = []
        for row in self.qubits:
            for q in row:
                if not isinstance(q, Data) or isinstance(q, SyData):
                    continue
                m_ls.append(Measurement(target=q, time=self.clock,xz='idk')) # TODO set xz parameter
        self.set_measurements(m_ls)

        string = ""
        for m in m_ls:
            string += f" {m.target.index}"
        all_other = self.all_other_qubits(string)
        return "M" + string + "\n" + f"DEPOLARIZE1({self.errors['measurement']})" + string + "\n" + f"DEPOLARIZE1({self.errors['idle']})" + all_other + "\n"

    def final_detectors(self):
        measurements = np.array(self.measurements)
        string = ""
        for h in self.hexagons:
            if(h.position[1] != 2):
                string += h.get_final_detectors(measurements,self.m_count)

        for m_ls in measurements[-1,:-1]: # bottom 3 body 
            if len(m_ls) == 0 or isinstance(m_ls[0].target, Data):
                continue

            m = m_ls[-1]
            x,y,_,_=m
            m1 = measurements[-1,x-1][-1]
            m2 = measurements[-1,x+1][-1]

            _,_,_,count1=m1
            _,_,_,count2=m2
            
            i = self.m_count - m.count
            i1 = self.m_count-count1
            i2 = self.m_count-count2

            string += f"DETECTOR rec[{-i}] rec[{-i1}] rec[{-i2}]\n"

        for m_ls in measurements[0]: # top 5 body
            if len(m_ls) == 0 or isinstance(m_ls[-1].target, SyData):
                continue
            
            m1 = m_ls[-1]
            x,y,t,count1=m1
            m2 = measurements[y,x-1][-1]
            m3 = measurements[y,x+1][-1]
            m4 = measurements[y+2,x-1][-1]
            m5 = measurements[y+2,x+1][-1]

            string += "DETECTOR"
            for m in [m1,m2,m3,m4,m5]:
                i = self.m_count - m.count
                string += f" rec[{-i}]"
            string += "\n"

        return string


    def observable(self):
        ms = []
        measurements = np.array(self.measurements)
        for m_ls in measurements[:,1]:
                if len(m_ls) > 0:
                    ms.append(m_ls[-1])
        string = "OBSERVABLE_INCLUDE(1)"
        for m in ms:
            i = self.m_count - m.count
            string += f" rec[{-i}]"
        return string

    def write_to_stim(self,):
        string = self.get_circuit()
        c = stim.Circuit(string)
        return c


    def __str__(self) -> str:

        string_grid = np.empty(shape=self.qubits.shape, dtype='<U3')

        for y, row in enumerate(self.qubits):
            for x, q in enumerate(row):
                if isinstance(q, Qubit):
                    string_grid[y,x] = f"{q.index:03}"
                else: 
                    string_grid[y,x] = "|-|"
        
        for h in self.hexagons:
            x, y = h.position[0], h.position[1]
            index = h.index
            string_grid[y,x] = f"H{index:02}"

        string = ""
        for row in string_grid:
            for s in row:
                string += f"{s} "
            string += "\n"
        
        return string