__author__ = 'Santiago Salas'
import cantera as ct
import csv
import numpy as np
import sys
from PyQt4 import QtGui
# define gui app and widget
app = QtGui.QApplication(sys.argv)
# Widget class with changed icon


class MainForm(QtGui.QWidget):

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.setGeometry(300, 100, 550, 500)
        self.setFixedWidth(400)
        self.setMinimumHeight(500)
        self.setWindowTitle('Icon')
        self.setWindowIcon(QtGui.QIcon('./utils/ch_pot.png'))

        self.table_properties = QtGui.QTableWidget(4, 2, self)
        self.table_state_funcs = QtGui.QTableWidget(6, 3, self)
        self.table_composition = QtGui.QTableWidget(1, 3, self)
        v_box = QtGui.QVBoxLayout(self)
        v_box.addStrut(5)
        v_box.addWidget(self.table_properties)
        v_box.addWidget(self.table_state_funcs)
        v_box.addWidget(self.table_composition)
        self.setLayout(v_box)

    def calculate_eq(self):
        ## CALCULATIONS
        # Instantiate
        gas = ct.Solution('gri30.cti')
        # Properties T,P,X
        gas.TPY = 300, 100000, 'CH4:1, O2:2'
        # Equilibrate @ constant TP
        print '\n\n Equilibrating... \n\n'
        gas.equilibrate('TP')
        print '\n\n ...Done... \n\n'
        return gas

    def print_eq(self, gas):
        ## PRINTING
        # Print to console
        print gas()
        # Save csv file in ./DATA/
        print '\n\n Saving csv... \n\n'
        csv_file = './DATA/equilibrium_simple_ex_01_result_ignorethis.csv'
        # Quickest output.
        with open(csv_file, 'w') as outfile:
            writer = csv.writer(outfile, lineterminator='\n')
            writer.writerow(['T (K)']+gas.species_names)
            writer.writerow([gas.T]+list(gas.X))
        # Redundant but more complete output.
        with open(csv_file, 'w') as outfile:
            sort_indexes = np.flipud(np.argsort(gas.Y))
            sorted_y = gas.Y[sort_indexes]
            sorted_x = gas.X[sort_indexes]
            sorted_names = np.array(gas.species_names)[sort_indexes]
            writer = csv.writer(outfile, lineterminator='\n')
            writer.writerow(['T', gas.T, 'K'])
            writer.writerow(['P', gas.P, 'Pa'])
            writer.writerow(['rho', gas.density, 'kg/m^3'])
            writer.writerow(['MW_mean', gas.mean_molecular_weight, 'amu'])
            writer.writerow([''])
            writer.writerow(['', '1 kg', '1 kmol', ''])
            writer.writerow(['H', gas.enthalpy_mass, gas.enthalpy_mole, 'J'])
            writer.writerow(['U', gas.int_energy_mass, gas.int_energy_mole, 'J'])
            writer.writerow(['S', gas.entropy_mass, gas.entropy_mole, 'J/K'])
            writer.writerow(['G', gas.gibbs_mass, gas.gibbs_mole, 'J'])
            writer.writerow(['c_p', gas.cp_mass, gas.cp_mole, 'J/K'])
            writer.writerow(['c_v', gas.cv_mass, gas.cv_mole, 'J/K'])
            writer.writerow([''])
            writer.writerow(['', 'X', 'Y', 'Chem. Pot. / RT'])

            self.table_composition.setRowCount(sort_indexes.size)
            for i in range(sort_indexes.size):
                writer.writerow([sorted_names[i]] + [sorted_x[i]] + [sorted_y[i]] +
                                [gas.chemical_potentials[i]/gas.T/ct.gas_constant])
                self.table_composition.setItem(i, 0, QtGui.QTableWidgetItem(str(sorted_x[i])))
                self.table_composition.setItem(i, 1, QtGui.QTableWidgetItem(str(sorted_y[i])))
                self.table_composition.setItem(i, 2,
                                               QtGui.QTableWidgetItem(
                                                   str(gas.chemical_potentials[i]/
                                                       gas.T/ct.gas_constant)))

            self.table_properties.setItem(0, 0, QtGui.QTableWidgetItem(str(gas.T)))
            self.table_properties.setItem(0, 1, QtGui.QTableWidgetItem('K'))
            self.table_properties.setItem(1, 0, QtGui.QTableWidgetItem(str(gas.P)))
            self.table_properties.setItem(1, 1, QtGui.QTableWidgetItem('Pa'))
            self.table_properties.setItem(2, 0, QtGui.QTableWidgetItem(str(gas.density)))
            self.table_properties.setItem(2, 1, QtGui.QTableWidgetItem('kg/m^3'))
            self.table_properties.setItem(3, 0, QtGui.QTableWidgetItem(
                str(gas.mean_molecular_weight)))
            self.table_properties.setItem(3, 1, QtGui.QTableWidgetItem('amu'))
            self.table_state_funcs.setItem(0, 0, QtGui.QTableWidgetItem(
                str(gas.enthalpy_mass)))
            self.table_state_funcs.setItem(1, 0, QtGui.QTableWidgetItem(
                str(gas.int_energy_mass)))
            self.table_state_funcs.setItem(2, 0, QtGui.QTableWidgetItem(
                str(gas.entropy_mass)))
            self.table_state_funcs.setItem(3, 0, QtGui.QTableWidgetItem(
                str(gas.gibbs_mass)))
            self.table_state_funcs.setItem(4, 0, QtGui.QTableWidgetItem(
                str(gas.cp_mass)))
            self.table_state_funcs.setItem(5, 0, QtGui.QTableWidgetItem(
                str(gas.cv_mass)))
            self.table_state_funcs.setItem(0, 1, QtGui.QTableWidgetItem(
                str(gas.enthalpy_mole)))
            self.table_state_funcs.setItem(1, 1, QtGui.QTableWidgetItem(
                str(gas.int_energy_mole)))
            self.table_state_funcs.setItem(2, 1, QtGui.QTableWidgetItem(
                str(gas.entropy_mole)))
            self.table_state_funcs.setItem(3, 1, QtGui.QTableWidgetItem(
                str(gas.gibbs_mole)))
            self.table_state_funcs.setItem(4, 1, QtGui.QTableWidgetItem(
                str(gas.cp_mole)))
            self.table_state_funcs.setItem(5, 1, QtGui.QTableWidgetItem(
                str(gas.cv_mole)))
            self.table_state_funcs.setItem(0, 2, QtGui.QTableWidgetItem('J'))
            self.table_state_funcs.setItem(1, 2, QtGui.QTableWidgetItem('J'))
            self.table_state_funcs.setItem(2, 2, QtGui.QTableWidgetItem('J/K'))
            self.table_state_funcs.setItem(3, 2, QtGui.QTableWidgetItem('J'))
            self.table_state_funcs.setItem(4, 2, QtGui.QTableWidgetItem('J/K'))
            self.table_state_funcs.setItem(5, 2, QtGui.QTableWidgetItem('J/K'))
            self.table_properties.setHorizontalHeaderLabels(['', ''])
            self.table_properties.setVerticalHeaderLabels(['T', 'P', 'rho', 'mean_MW'])
            self.table_state_funcs.setHorizontalHeaderLabels(['1 kg', '1 kmol', ''])
            self.table_state_funcs.setVerticalHeaderLabels(['H', 'U', 'S', 'G', 'Cp', 'Cv'])
            self.table_composition.setHorizontalHeaderLabels(['X', 'Y', 'Chem. Pot. / RT'])
            self.table_composition.setVerticalHeaderLabels(sorted_names)
        print '\n\n ...Done... \n\n'
icon = MainForm()
icon.print_eq(icon.calculate_eq())
icon.show()
icon.setWindowTitle('CH4 / O2 Equilaibrium composition')
# show the main widget, enter 'main loop'

# exit 'main loop'
sys.exit(app.exec_())