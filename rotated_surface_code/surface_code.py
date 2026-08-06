import stim
import numpy as np

class Qubit():
    def __init__(self, index: int, position: tuple):
        self.index = index
        self.pos = np.array([position[0], position[1]])

    def __str__(self):
        x, y = self.pos
        i = self.index
        return f"QUBIT_COORDS({x},{y}) {i}"
    
    def get_pos(self):
        return (self.pos[0], self.pos[1])
    
    @staticmethod
    def find_at(coord, qubit_ls, reverse=False):
        for q in qubit_ls:
            if(q.get_pos()[0] == coord[0] and q.get_pos()[1] == coord[1]):
                return q
        return None
    
    def where_in(self, qubit_ls: list, reverse=False):
        indices = []
        for i, q in enumerate(qubit_ls):
            if q == self:
                if reverse:
                    indices.append(-len(qubit_ls) + i)
                else:
                    indices.append(i)
        return indices

    def __eq__(self, other_qubit):
        if(self.index == other_qubit.index) and (self.pos[0] == other_qubit.pos[0] and self.pos[1] == other_qubit.pos[1]):
            return True
        return False
        
class Zancilla(Qubit):

    schedule_coords = np.array([[0.5,0.5] , [0.5,-0.5] , [-0.5,0.5] , [-0.5,-0.5]])

    def __init__(self, index:int, position:tuple):
        super().__init__(index, position)
        
    def set_schedule(self, data_qubits: list[Qubit]):
        schedule = []
        coords = self.schedule_coords + np.array( [self.pos for i in range(4)] )
        for coord in coords:
            schedule.append(Qubit.find_at(coord, data_qubits))
        self.schedule = schedule
    
class Xancilla(Qubit):

    schedule_coords = np.array([[0.5,0.5] , [-0.5,0.5] , [0.5,-0.5] , [-0.5,-0.5]])

    def __init__(self, index:int, position:tuple):
        super().__init__(index, position)
        self.schedule = None

    def set_schedule(self, data_qubits: list[Qubit]):
        schedule = []
        coords = self.schedule_coords + np.array( [self.pos for i in range(4)] )
        for coord in coords:
            schedule.append(Qubit.find_at(coord, data_qubits))
        self.schedule = schedule


class GenerateSurfaceCode():

    def __init__(self, distance, reps=24, filename="surface_code.stim", noise=0.01):
        self.string = ""
        self.dist = distance
        self.filename = filename
        self.reps = reps

        self.single_error = noise/10
        self.double_error = noise
        self.idle_error = noise*2
        self.reset_error = noise*2
        self.measure_error = noise*5
        
        self.data_qubits = self.initialise_data_qubits()
        self.z_ancillas= self.initialise_z_ancillas()
        self.x_ancillas = self.initialise_x_ancillas()

        self.set_schedules() 
        self.measurements = []

        self.write_circuit()
         
    def initialise_data_qubits(self):
        """
        Initialises nxn grid of data qubits with positions starting at (1,1) and on integer coords
        """
        ls = []
        for i in range(self.dist**2):
            x = i%self.dist + 1
            y = i//self.dist+1
            ls.append( Qubit(index=i, position=(x,y)) )
        return ls

    def initialise_z_ancillas(self):
        """
        Initialises placement of z ancillas with position starting at (2.5, 0.5) and going left to right - up to down
        """
        nd = self.dist**2 # number of data qubits (for indexing)
        za_ls = []        
        
        for y in range(1, self.dist):
            for x in range(0, self.dist+1):
                x_pos = x+0.5
                y_pos = y+0.5
                if( (x+y)%2==0 ):
                    za_ls.append(Zancilla( nd+len(za_ls), (x_pos,y_pos) ))
        return za_ls
        
    def initialise_x_ancillas(self):
        """
        Initialises placement of x ancillas with position starting at (0.1, 1.5) and going left to right - up to down
        """
        nd = self.dist**2 # number of data qubits (for indexing)
        nz = len(self.z_ancillas)
        xa_ls = []
    
        for y in range(0, self.dist+1):
            for x in range(1, self.dist):
                x_pos = x+0.5
                y_pos = y+0.5
                if( (x+y)%2 == 1 ):
                    xa_ls.append(Xancilla( index=nd+nz+len(xa_ls), position=(x_pos,y_pos) ))
        return xa_ls
        
    def set_schedules(self):
        for za in self.z_ancillas:
            za.set_schedule(self.data_qubits)
        for xa in self.x_ancillas:
            xa.set_schedule(self.data_qubits)

    def write_qubits(self) -> str:
        """
        writes qubits into string with proper stim syntax to intitialise circuit
        """
        string = ""
        qubits = self.data_qubits + self.z_ancillas + self.x_ancillas
        for q in qubits:
            string += str(q) + "\n"
        self.string += string + "\n"
        return string
    
    def reset_all(self):
        string = "R "
        qubits = self.data_qubits + self.z_ancillas + self.x_ancillas
        for q in qubits:
            string += f" {q.index}"
        string += self.write_reset_error(qubits)
        self.string += string + "\n"
        self.tick()
    
    def write_hadamards(self):
        string = "H "
        qubits = self.x_ancillas
        for q in qubits:
            string += f" {q.index}"
        string += self.write_depolarise_error(qubits)
        self.string += string + "\n"
        self.tick()
    
    def write_CNOTs(self):
        string = ""

        for i in range(4):
            qubits = []
            string += "CX"

            for z in self.z_ancillas:
                if isinstance(z.schedule[i], Qubit):
                    control = z.schedule[i]
                    test = z
                    string += f" {control.index} {test.index}"

                    qubits.append(control)
                    qubits.append(test)

            for x in self.x_ancillas:
                if isinstance(x.schedule[i], Qubit):
                    test = x.schedule[i]
                    control = x
                    string += f" {control.index} {test.index}"

                    qubits.append(control)
                    qubits.append(test)
            
            string += self.write_depolarise_error(qubits, n=2)
            string += "\nTICK\n"

        self.string += string + "\n"

    def tick(self):
        self.string += "\nTICK\n"
    
    def measure_ancillas(self, n=0):
        
        qubits = self.z_ancillas + self.x_ancillas
        string = ""
        string += self.write_measure_error(qubits)

        string += "\nM"
        for q in qubits:
            string += f" {q.index}"
            self.measurements.append(q)

        if n == 0:
            string += "\nTICK\nR"

            for q in qubits:
                string += f" {q.index}"
        
            string += "\n" + self.write_reset_error(qubits)

        self.string += string + "\n"
        if n==0:
            self.tick()

    def write_detectors(self, condition=1):
        string = ""
        if(condition == 0):
            for z in self.z_ancillas: # first run
                x, y = z.get_pos()
                index = z.where_in(self.measurements, reverse=True)[-1]
                string += f"DETECTOR({x}, {y}) rec[{index}]\n" 

        elif(condition == 1): # in repeat block
            for a in self.z_ancillas + self.x_ancillas:
                x,y = a.get_pos()
                indices = a.where_in(self.measurements, reverse=True)
                index1, index2 = indices[-1], indices[-2]
                string += f"DETECTOR({x}, {y}) rec[{index1}] rec[{index2}]\n" 

        elif condition == 2: # last run
            for z in self.z_ancillas:
                x,y=z.get_pos()
                string += f"DETECTOR({x}, {y})"
                for d in z.schedule:
                    if isinstance(d, Qubit):
                        index = d.where_in(self.measurements, reverse=True)[-1]
                        string += f" rec[{index}]"
                index = z.where_in(self.measurements, reverse=True)[-1]
                string += f" rec[{index}]"
                string += "\n"

        self.string += string
    
    def open_repeat(self):
        string = f"REPEAT {self.reps} {{\n" 
        self.string += string 
        
    def close_repeat(self):
        self.string += f"}}\n"

    def write_reset_error(self, qubits):
        string = f"\nX_ERROR({self.reset_error})"
        for q in qubits:
            string += f" {q.index}"

        string += f"\nDEPOLARIZE1({self.idle_error})"
        other_qubits = self.data_qubits + self.x_ancillas + self.z_ancillas
        for q_other in other_qubits:
            if q_other not in qubits:
                string += f" {q_other.index}"

        return string + "\n"

    def write_measure_error(self, qubits):
        string = f"\nX_ERROR({self.measure_error})"
        for q in qubits:
            string += f" {q.index}"

        string += f"\nDEPOLARIZE1({self.idle_error})"
        other_qubits = self.data_qubits + self.x_ancillas + self.z_ancillas
        for q_other in other_qubits:
            if q_other not in qubits:
                string += f" {q_other.index}"

        return string + "\n"

    def write_depolarise_error(self, qubits, n=1):
        if(n == 1):
            string = f"\n DEPOLARIZE{n}({self.single_error})"
        else:
            string = f"\n DEPOLARIZE{n}({self.double_error})"

        for q in qubits:
            string += f" {q.index}"
        return string + "\n"

    def measure_data(self):
        string = ""
        string += self.write_measure_error(self.data_qubits)

        string += "\nM"
        for d in self.data_qubits:
            self.measurements.append(d)
            string += f" {d.index}"
        self.string += string + "\n"

    def write_observable(self):
        string = ""
        string += f"\nOBSERVABLE_INCLUDE(0)" # top row
        indices = [i for i in range(0,self.dist)]
        for d in np.array(self.data_qubits)[indices]:
            index = d.where_in(self.measurements, reverse=True)[-1]
            string += f" rec[{index}]"

        self.string += string

    def write_circuit(self):
        self.write_qubits()       
        self.reset_all()          
        self.write_hadamards()    
        self.write_CNOTs()        
        self.write_hadamards() 
        self.measure_ancillas()
        self.write_detectors(0)
        
        if self.reps > 0:
            self.open_repeat()

            self.write_hadamards()
            self.write_CNOTs()  
            self.write_hadamards()
            self.measure_ancillas()
            self.write_detectors()

            self.close_repeat()

        self.write_hadamards()
        self.write_CNOTs()  
        self.write_hadamards()
        self.measure_ancillas(n=1)
        self.write_detectors()

        self.measure_data()
        self.write_detectors(2)
        self.write_observable()
        
    def get_circuit(self):
        return stim.Circuit(self.string)
