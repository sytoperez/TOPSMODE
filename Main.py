#==============================================================================
# author          : Alfonso Pérez-Garrido
# date            : 01-10-2023
# version         : 0.1
# python_version  : 3.7
#==============================================================================


import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from PIL import Image, ImageTk
import os
import sys
import csv
import TOPS
import TOPS_Ctbr



def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


class DescriptorsFrame(ttk.Labelframe):

    def __init__(self, parent, name, column, row, full_version=True):

        ttk.Labelframe.__init__(self, parent, text=name)
        self.grid(column=column, row=row, columnspan=10, sticky=(tk.W, tk.E), padx=5, pady=5)

        self._full_version = full_version

        self.des_type = tk.StringVar(value='bond')

        self.des_std = tk.BooleanVar(value=False)
        self.des_dip = tk.BooleanVar(value=False)
        self.des_dip2 = tk.BooleanVar(value=False)
        self.des_hyd = tk.BooleanVar(value=False)
        self.des_pls = tk.BooleanVar(value=False)
        self.des_mr = tk.BooleanVar(value=False)
        self.des_pol = tk.BooleanVar(value=False)
        self.des_vdw = tk.BooleanVar(value=False)
        self.des_chg = tk.BooleanVar(value=False)
        self.des_atw = tk.BooleanVar(value=False)
        self.des_ar2 = tk.BooleanVar(value=False)
        self.des_api2 = tk.BooleanVar(value=False)
        self.des_a2h = tk.BooleanVar(value=False)
        self.des_b2h = tk.BooleanVar(value=False)
        self.des_b2o = tk.BooleanVar(value=False)
        self.des_l16 = tk.BooleanVar(value=False)

        self.des_std_b = tk.BooleanVar(value=True)
        self.des_dip_b = tk.BooleanVar(value=True)
        self.des_dip2_b = tk.BooleanVar(value=True)
        self.des_hyd_b = tk.BooleanVar(value=True)
        self.des_pls_b = tk.BooleanVar(value=True)
        self.des_mr_b = tk.BooleanVar(value=True)
        self.des_pol_b = tk.BooleanVar(value=True)
        self.des_vdw_b = tk.BooleanVar(value=True)
        self.des_chg_b = tk.BooleanVar(value=True)
        self.des_atw_b = tk.BooleanVar(value=True)
        self.des_ar2_b = tk.BooleanVar(value=True)
        self.des_api2_b = tk.BooleanVar(value=True)
        self.des_a2h_b = tk.BooleanVar(value=True)
        self.des_b2h_b = tk.BooleanVar(value=True)
        self.des_b2o_b = tk.BooleanVar(value=True)
        self.des_l16_b = tk.BooleanVar(value=True)
        self._value = tk.StringVar(value=str(15))
        ttk.Radiobutton(self, text='Atomic contribution', name='type_ato', value='ato', variable=self.des_type, command=self.radio_atom).grid(
            column=0, row=0, sticky=(tk.W), padx=5, pady=5)
        ttk.Radiobutton(self, text='Bond contributions', name='type_bond', value='bond', variable=self.des_type, command=self.radio_bond).grid(
            column=6, row=0, sticky=(tk.E,tk.W), padx=5, pady=5)

        if self._full_version:
            ttk.Checkbutton(self, variable=self.des_std, name='chk_std_a', text='Standard distance',
                            command=self.ato_des_checked).grid(column=0, row=3, sticky=(tk.W), padx=5,
                                                               pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_dip, name='chk_dip_a', text='Dipole moment',
                            command=self.ato_des_checked).grid(column=0, row=4, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_dip2, name='chk_dip2_a', text='Dipole moment 2',
                            command=self.ato_des_checked).grid(column=0, row=5, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_hyd, name='chk_hyd_a', text='Hydrophobicity',
                            command=self.ato_des_checked).grid(column=0, row=6, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_pls, name='chk_pls_a', text='Polar Surface',
                            command=self.ato_des_checked).grid(column=0, row=7, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_mr, name='chk_mr_a', text='Molar refractivity',
                            command=self.ato_des_checked).grid(column=0, row=8, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_pol, name='chk_pol_a', text='Polarizability',
                            command=self.ato_des_checked).grid(column=0, row=9, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_vdw, name='chk_vdw_a', text='Van der Waal radius',
                            command=self.ato_des_checked).grid(column=0, row=10, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_chg, name='chk_chg_a', text='Gasteiger chargue',
                            command=self.ato_des_checked).grid(column=1, row=3, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_atw, name='chk_atw_a', text='Atomic weight',
                            command=self.ato_des_checked).grid(column=1, row=4, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_ar2, name='chk_ar2_a', text='Abraham-R2',
                            command=self.ato_des_checked).grid(column=1, row=5, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_api2, name='chk_api2_a', text='Abraham-pi2H',
                            command=self.ato_des_checked).grid(column=1, row=6, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_a2h, name='chk_a2h_a', text='Abraham-a2H',
                            command=self.ato_des_checked).grid(column=1, row=7, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_b2h, name='chk_b2h_a', text='Abraham-B2H',
                            command=self.ato_des_checked).grid(column=1, row=8, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_b2o, name='chk_b2o_a', text='Abraham-B2O',
                            command=self.ato_des_checked).grid(column=1, row=9, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_l16, name='chk_l16_a', text='Abraham-L16',
                            command=self.ato_des_checked).grid(column=1, row=10, sticky=(tk.W), padx=5,
                                                                 pady=(1, 0))

            ttk.Checkbutton(self, variable=self.des_std_b, name='chk_std_b', text='Standard distance',
                            command=self.bond_des_checked).grid(column=6, row=3, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_dip_b, name='chk_dip_b', text='Dipole moment',
                            command=self.bond_des_checked).grid(column=6, row=4, sticky=( tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_dip2_b, name='chk_dip2_b', text='Dipole moment 2',
                            command=self.bond_des_checked).grid(column=6, row=5, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_hyd_b, name='chk_hyd_b', text='Hydrophobicity',
                            command=self.bond_des_checked).grid(column=6, row=6, sticky=( tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_pls_b, name='chk_pls_b', text='Polar Surface',
                            command=self.bond_des_checked).grid(column=6, row=7, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_mr_b, name='chk_mr_b', text='Molar refractivity',
                            command=self.bond_des_checked).grid(column=6, row=8, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_pol_b, name='chk_pol_b', text='Polarizability',
                            command=self.bond_des_checked).grid(column=6, row=9, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_vdw_b, name='chk_vdw_b', text='Van der Waal radius',
                            command=self.bond_des_checked).grid(column=6, row=10, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_chg_b, name='chk_chg_b', text='Gasteiger chargue',
                            command=self.bond_des_checked).grid(column=7, row=3, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_atw_b, name='chk_atw_b', text='Atomic weight',
                            command=self.bond_des_checked).grid(column=7, row=4, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_ar2_b, name='chk_ar2_b', text='Abraham-R2',
                            command=self.bond_des_checked).grid(column=7, row=5, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_api2_b, name='chk_api2_b', text='Abraham-pi2H',
                            command=self.bond_des_checked).grid(column=7, row=6, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_a2h_b, name='chk_a2h_b', text='Abraham-a2H',
                            command=self.bond_des_checked).grid(column=7, row=7, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_b2h_b, name='chk_b2h_b', text='Abraham-B2H',
                            command=self.bond_des_checked).grid(column=7, row=8, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_b2o_b, name='chk_b2o_b', text='Abraham-B2O',
                            command=self.bond_des_checked).grid(column=7, row=9, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))
            ttk.Checkbutton(self, variable=self.des_l16_b, name='chk_l16_b', text='Abraham-L16',
                            command=self.bond_des_checked).grid(column=7, row=10, sticky=(tk.E,tk.W), padx=5,
                                                                 pady=(1, 0))

        self.columnconfigure(0, pad=80)
    def radio_atom(self):
        self.des_std.set(value=True)
        self.des_dip.set(value=True)
        self.des_dip2.set(value=True)
        self.des_hyd.set(value=True)
        self.des_pls.set(value=True)
        self.des_mr.set(value=True)
        self.des_pol.set(value=True)
        self.des_vdw.set(value=True)
        self.des_chg.set(value=True)
        self.des_atw.set(value=True)
        self.des_ar2.set(value=True)
        self.des_api2.set(value=True)
        self.des_a2h.set(value=True)
        self.des_b2h.set(value=True)
        self.des_b2o.set(value=True)
        self.des_l16.set(value=True)
        self.des_std_b.set(value=False)
        self.des_dip_b.set(value=False)
        self.des_dip2_b.set(value=False)
        self.des_hyd_b.set(value=False)
        self.des_pls_b.set(value=False)
        self.des_mr_b.set(value=False)
        self.des_pol_b.set(value=False)
        self.des_vdw_b.set(value=False)
        self.des_chg_b.set(value=False)
        self.des_atw_b.set(value=False)
        self.des_ar2_b.set(value=False)
        self.des_api2_b.set(value=False)
        self.des_a2h_b.set(value=False)
        self.des_b2h_b.set(value=False)
        self.des_b2o_b.set(value=False)
        self.des_l16_b.set(value=False)

    def radio_bond(self):
        self.des_std.set(value=False)
        self.des_dip.set(value=False)
        self.des_dip2.set(value=False)
        self.des_hyd.set(value=False)
        self.des_pls.set(value=False)
        self.des_mr.set(value=False)
        self.des_pol.set(value=False)
        self.des_vdw.set(value=False)
        self.des_chg.set(value=False)
        self.des_atw.set(value=False)
        self.des_ar2.set(value=False)
        self.des_api2.set(value=False)
        self.des_a2h.set(value=False)
        self.des_b2h.set(value=False)
        self.des_b2o.set(value=False)
        self.des_l16.set(value=False)
        self.des_std_b.set(value=True)
        self.des_dip_b.set(value=True)
        self.des_dip2_b.set(value=True)
        self.des_hyd_b.set(value=True)
        self.des_pls_b.set(value=True)
        self.des_mr_b.set(value=True)
        self.des_pol_b.set(value=True)
        self.des_vdw_b.set(value=True)
        self.des_chg_b.set(value=True)
        self.des_atw_b.set(value=True)
        self.des_ar2_b.set(value=True)
        self.des_api2_b.set(value=True)
        self.des_a2h_b.set(value=True)
        self.des_b2h_b.set(value=True)
        self.des_b2o_b.set(value=True)
        self.des_l16_b.set(value=True)

    def bond_des_checked(self):
        self.des_type.set(value='bond')
        self.des_std.set(value=False)
        self.des_dip.set(value=False)
        self.des_dip2.set(value=False)
        self.des_hyd.set(value=False)
        self.des_pls.set(value=False)
        self.des_mr.set(value=False)
        self.des_pol.set(value=False)
        self.des_vdw.set(value=False)
        self.des_chg.set(value=False)
        self.des_atw.set(value=False)
        self.des_ar2.set(value=False)
        self.des_api2.set(value=False)
        self.des_a2h.set(value=False)
        self.des_b2h.set(value=False)
        self.des_b2o.set(value=False)
        self.des_l16.set(value=False)

    def ato_des_checked(self):
        self.des_type.set(value='ato')
        self.des_std_b.set(value=False)
        self.des_dip_b.set(value=False)
        self.des_dip2_b.set(value=False)
        self.des_hyd_b.set(value=False)
        self.des_pls_b.set(value=False)
        self.des_mr_b.set(value=False)
        self.des_pol_b.set(value=False)
        self.des_vdw_b.set(value=False)
        self.des_chg_b.set(value=False)
        self.des_atw_b.set(value=False)
        self.des_ar2_b.set(value=False)
        self.des_api2_b.set(value=False)
        self.des_a2h_b.set(value=False)
        self.des_b2h_b.set(value=False)
        self.des_b2o_b.set(value=False)
        self.des_l16_b.set(value=False)

    def get_selected_descriptors(self):

        output = []
        if self._full_version:
            if self.des_type.get() == 'bond':
                if self.des_std_b.get(): output.append('Std')
                if self.des_dip_b.get(): output.append('Dip')
                if self.des_dip2_b.get(): output.append('Dip2')
                if self.des_hyd_b.get(): output.append('Hyd')
                if self.des_pls_b.get(): output.append('Pols')
                if self.des_mr_b.get(): output.append('Mol')
                if self.des_pol_b.get(): output.append('Pol')
                if self.des_vdw_b.get(): output.append('Van')
                if self.des_chg_b.get(): output.append('Gas')
                if self.des_atw_b.get(): output.append('Ato')
                if self.des_ar2_b.get(): output.append('Ab-R2')
                if self.des_api2_b.get(): output.append('Ab-pi2H')
                if self.des_a2h_b.get(): output.append('Ab-sumA2H')
                if self.des_b2h_b.get(): output.append('Ab-sumB2H')
                if self.des_b2o_b.get(): output.append('Ab-sumB20')
                if self.des_l16_b.get(): output.append('Ab-logL16')
            else:
                if self.des_std.get(): output.append('Std')
                if self.des_dip.get(): output.append('Dip')
                if self.des_dip2.get(): output.append('Dip2')
                if self.des_hyd.get(): output.append('Hyd')
                if self.des_pls.get(): output.append('Pols')
                if self.des_mr.get(): output.append('Mol')
                if self.des_pol.get(): output.append('Pol')
                if self.des_vdw.get(): output.append('Van')
                if self.des_chg.get(): output.append('Gas')
                if self.des_atw.get(): output.append('Ato')
                if self.des_ar2.get(): output.append('Ab-R2')
                if self.des_api2.get(): output.append('Ab-pi2H')
                if self.des_a2h.get(): output.append('Ab-sumA2H')
                if self.des_b2h.get(): output.append('Ab-sumB2H')
                if self.des_b2o.get(): output.append('Ab-sumB20')
                if self.des_l16.get(): output.append('Ab-logL16')
        else:
            output = ['Std', 'Dip', 'Dip2', 'Hyd', 'Pols', 'Mol', 'Pol', 'Van', 'Gas', 'Ato', 'Ab-R2',
                             'Ab-pi2H', 'Ab-sumA2H', 'Ab-sumB2H', 'Ab-sumB20', 'Ab-logL16']

        return output

    def get_selected_descriptors_with_type(self):
        descriptors = self.get_selected_descriptors()
        return [d + '_' + self.des_type.get() for d in descriptors]


class NumberFrame(ttk.Labelframe):

    def __init__(self, parent, column, row, name="Max Order of descriptors:"):

        ttk.Labelframe.__init__(self, parent, text=name)
        self.grid(column=column, row=row, columnspan=10, sticky=(tk.W, tk.E), padx=5, pady=5)

        self._value = tk.StringVar(value=str(15))
        tk.Spinbox(self, from_=1, to=15, textvariable=self._value, width=5).grid(column=0, row=2, columnspan=10, sticky=(tk.W, tk.E), padx=5, pady=(1, 0))

    def get_value(self):
        return int(self._value.get())


class Sdfreader(ttk.Labelframe):
    def __select_sdf_path(self):
        self.sdf_path.set(filedialog.askopenfilename(filetypes=[('SDF files', '*.sdf')]))

    def __select_property_file_path(self):
        self.property_file_path.set(filedialog.askopenfilename(filetypes=[('Property text file', '*.txt')],
                                                               initialdir=os.path.dirname(self.sdf_path.get())))

    def __remove_forbidden_y(self, filename):
        # remove non-numeric items from y file (with compounds properties)

        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        lines = open(filename).readlines()
        # add header to the output
        output = [lines[0]]
        for line in lines[1:]:
            if is_number(line.strip().split('\t')[1]):
                output.append(line)
        open(filename, 'wt').writelines(output)

    def __add_sdf_path_to_history(self, sdf_path):
        if sdf_path is None or sdf_path == "":
            return None
        max_lines = 10
        new_line = sdf_path + '\n'
        hist = os.path.join(get_script_path(), 'history.txt')
        if os.path.isfile(hist):
            lines = open(hist).readlines()
            # remove newly added line if it is present in lines
            if new_line in lines:
                lines.remove(new_line)
            # keep only allowed number of lines - 1
            if len(lines) >= max_lines:
                lines = lines[-(max_lines-1):]
            lines.append(new_line)
            open(hist, 'wt').writelines(lines)
        else:
            open(hist, 'wt').write(new_line)

    def __read_sdf_history(self):
        file_name = os.path.join(get_script_path(), 'history.txt')
        if os.path.isfile(file_name):
            lines = open(file_name).readlines()
            lines = [line.strip() for line in lines]
            lines.reverse()
            return lines

    def __sdf_path_changed(self, varname, elementname, mode):
        field_names = self.__read_sdf_field_names(self.sdf_path.get())
        if field_names:
            self.children['sdf_label_frame'].children['activity_field_name'].configure(values=field_names)
            self.children['sdf_label_frame'].children['activity_field_name'].set(value=field_names[0])
            #self.children['optional_label_frame'].children['inner_frame'].children['sdf_id_field_name'].configure(
            #    values=field_names)
            #self.children['optional_label_frame'].children['inner_frame'].children['sdf_id_field_name'].set(
            #    value=field_names[0])
        self.__add_sdf_path_to_history(self.sdf_path.get())
        self.children['sdf_label_frame'].children['sdf_path_combobox'].configure(values=self.__read_sdf_history())

        # update list of models to plot
        #self.master.children['tab_3']._show_models_list()

    def __read_sdf_field_names(self, fname):
        field_names = []
        if not os.path.isfile(fname):
            messagebox.showerror('ERROR!', "Specified file name doesn't exist.")
            return field_names
        with open(fname) as f:
            line = f.readline().rstrip()
            while line != '$$$$':
                # one or two spaces between > and < can be possible
                if line.startswith('>  <') and line.endswith('>'):
                    field_names.append(line[4:-1])
                if line.startswith('> <') and line.endswith('>'):
                    field_names.append(line[3:-1])
                line = f.readline().rstrip()
        return field_names

    def __init__(self, parent, name, column, row, full_version=True):
        ttk.Labelframe.__init__(self, parent, text=name, name='mol_input')
        self.grid(column=column, row=row, columnspan=10, sticky=(tk.W, tk.E), padx=5, pady=5)
        self._full_version = full_version
        self.sdf_path = tk.StringVar()
        self.sdf_id_field_name = tk.StringVar()
        self.compound_names = tk.StringVar(value='gen')
        self.activity_field_name = tk.StringVar()

        frame = ttk.Labelframe(self, text='SDF with compounds', name='sdf_label_frame')
        frame.grid(column=0, row=2, sticky=(tk.E, tk.W), columnspan=10, padx=5, pady=5)
        ttk.Label(frame, text='Path to SDF-file').grid(column=0, row=2, sticky=(tk.W, tk.S), padx=5)
        ttk.Label(frame, text='activity field name').grid(column=2, row=2, sticky=(tk.W), padx=5)
        ttk.Combobox(frame, name='sdf_path_combobox', width=70, textvariable=self.sdf_path,
                     values=self.__read_sdf_history()).grid(column=0, row=3, sticky=(tk.W, tk.E), padx=5, pady=(0, 5))
        ttk.Button(frame, text='Browse...', command=self.__select_sdf_path).grid(column=1, row=3, sticky=(tk.W), padx=5,
                                                                                 pady=(0, 5))
        ttk.Combobox(frame, name='activity_field_name', width=20, textvariable=self.activity_field_name, state='readonly').grid(column=2, row=3, sticky=(tk.W), padx=5, pady=(0, 5))

        frame = ttk.Labelframe(self, text='Optional. Compound names. External text file with compound property values',
                               name='optional_label_frame')
        frame.grid(column=10, row=2, sticky=(tk.E, tk.W), columnspan=5, padx=5, pady=5)

        ttk.Radiobutton(frame, text='Automatically generate compound names', variable=self.compound_names,
                        value='gen').grid(column=0, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Use compound titles from SDF file', variable=self.compound_names,
                        value='title').grid(column=0, row=1, sticky=(tk.W), padx=5, pady=1)

        self.sdf_path.trace('w', self.__sdf_path_changed)
        #self.property_field_name.trace('w', self.__act_field_name_changed)

class Tab_1(ttk.Frame):

    def __build_descriptors(self):

        if self.path.get() == '':
            messagebox.showerror('ERROR!', 'Specify sdf filename and path.')
            return

        print("Descriptors calculation started. Please wait it can take some time")

        TOPS.main_params(in_fname=self.path.get(), out_fname='TOPSMODE.csv',
                         des_names=self.des_frame.get_selected_descriptors(),
                         activity_field_set=self.activity.get(),
                         des_type=self.des_frame.des_type.get(), numero=self.number_count.get_value(), id_field_set=self.nombres.get(), verbose=0)

        print("Descriptors calculation finished")

    def __init__(self, parent, tab_name, path, nombres, act_field):

        ttk.Frame.__init__(self, parent, name='tab_1')
        self.path = path
        self.nombres = nombres
        self.des_frame = DescriptorsFrame(self, 'Descriptors', 0, 1, True)
        self.number_count = NumberFrame(self, 0, 0)
        self.activity = act_field
        ttk.Button(self, text='Calculate descriptors', command=self.__build_descriptors).grid(column=0, row=2,
                                                                                              columnspan=5,
                                                                                              sticky=(tk.W, tk.E))

        self.columnconfigure(0, weight=1)

        parent.add(self, text=tab_name)


class Tab_2(ttk.Frame):
    def __select_csv_path(self):
        self.csv_path.set(filedialog.askopenfilename(filetypes=[('CSV files', '*.csv')]))

    def __read_csv_history(self):
        file_name = os.path.join(get_script_path(), 'history_csv.txt')
        if os.path.isfile(file_name):
            lines = open(file_name).readlines()
            lines = [line.strip() for line in lines]
            lines.reverse()
            return lines

    def __contr_changed(self, varname, elementname, mode):
        if (self.tipo=='ato'):
            camino = os.path.join(os.path.dirname(self.path.get()),
                                  "TOPSMODE/Contr_ato" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "_" + self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "\\modelo_" + str(self.n_mod) + "\\" + str(
                                      self.contr_name.get()) + "_" + str(
                                      self.n_mol) + '.png')
        else:
            camino = os.path.join(os.path.dirname(self.path.get()),
                                  "TOPSMODE/Contr_" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "_" + self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "\\modelo_" + str(self.n_mod) + "\\" + str(
                                      self.contr_name.get()) + "_" + str(
                                      self.n_mol) + '.png')

        if os.path.exists(camino):
            self.original = Image.open(camino)
        else:
            self.original = Image.open(os.path.join(get_script_path(), 'No_mol.png'))
        resized = self.original.resize((400, 400), Image.ANTIALIAS)
        self.image = ImageTk.PhotoImage(resized)  # Keep a reference, prevent GC
        self.l2.config(image=self.image)
        self.label.config(text=self.names[self.n_mol - 1])

    def __csv_path_changed(self,varname, elementname, mode):
        if self.csv_path.get() == '':
            messagebox.showerror('ERROR!', 'Specify sdf filename and path.')
            return
        self.__add_csv_path_to_history(self.csv_path.get())
        self.children['csv_label_frame'].children['csv_path_combobox'].configure(values=self.__read_csv_history())
        i = 0
        self.modelos=[]
        self.tree.delete(self.tree.get_children())
        with open(self.csv_path.get()) as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                i = i + 1
                self.tree.insert('', 'end', values=('modelo_'+str(i), row['n'], row['variables'], row['coeff'], "--"),
                           tags='checked')
                self.modelos.append(['modelo_'+str(i), int(row['n']), str(row['variables']), str(row['coeff'])])
        print('Contribution calculation started. Please wait it can take some time')
        v = self.modelos[0][2].split('|')
        v[0] = v[0].replace('(', '|')
        v[0] = v[0].replace(')', '|')
        clase = v[0].split('|')
        if clase[0] == 'u':
            self.tipo ='bond'
        else:
            self.tipo ='ato'
        print(self.tipo)
        self.names = TOPS_Ctbr.main_params(in_fname=self.path.get(), in_model=self.csv_path.get(),
                                           out_fname='TOPSMODE_contr', id_field_set=self.nombres.get(),
                                           type_set=self.tipo, linear_set='lin', data_only=self.datos.get(), verbose=0)
        self.tree.insert('', 'end', values=('modelo_' + str(i+1), 'consenso', '-', '-', "--"), tags='checked')
        self.modelos.append(['modelo_' + str(i+1), 'consenso', '-', '-'])
        print('Calculation finished')


    def __add_csv_path_to_history(self, csv_path):
        if csv_path is None or csv_path == "":
            return None
        max_lines = 10
        new_line = csv_path + '\n'
        hist = os.path.join(get_script_path(), 'history_csv.txt')
        if os.path.isfile(hist):
            lines = open(hist).readlines()
            # remove newly added line if it is present in lines
            if new_line in lines:
                lines.remove(new_line)
            # keep only allowed number of lines - 1
            if len(lines) >= max_lines:
                lines = lines[-(max_lines-1):]
            lines.append(new_line)
            open(hist, 'wt').writelines(lines)
        else:
            open(hist, 'wt').write(new_line)

    def __select_save_file_path(self):
        self.save_file_path.set(filedialog.asksaveasfilename(filetypes=[('PNG files', '*.png')],
                                                             defaultextension='.png',
                                                         initialdir=os.path.dirname(self.path.get())))

    def __init__(self, parent, tab_name, path, nombres):

        ttk.Frame.__init__(self, parent, name='tab_2')
        self.path = path
        self.nombres = nombres
        self.csv_path = tk.StringVar()
        self.tipo = 'bond'
        self.modelos = []
        self.contr_overall = tk.BooleanVar(value=True)
        self.modelo_check = False
        self.n_mol = 1
        self.n_mod = 0
        self.names = []
        self.datos = tk.StringVar(value='total')
        self.contr_name = tk.StringVar(value='total')
        frame = ttk.Labelframe(self, text='Select models list to visualise', name='csv_label_frame')
        frame.grid(column=0, row=0, columnspan=10, sticky=(tk.E, tk.W), padx=5, pady=5)
        # select models
        ttk.Label(frame, text='Path to CSV-file').grid(column=0, row=0, sticky=(tk.W, tk.S), padx=5)
        ttk.Combobox(frame, name='csv_path_combobox', width=70, textvariable=self.csv_path,
                     values=self.__read_csv_history()).grid(column=0, row=1, sticky=(tk.W, tk.E), padx=5, pady=(0, 5))
        ttk.Button(frame, text='Browse...', command=self.__select_csv_path).grid(column=1, row=1, sticky=(tk.W), padx=5,
                                                                                 pady=(0, 5))
        ttk.Radiobutton(frame, text='Only Data', variable=self.datos,
                        value='data').grid(column=2, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Data+Structure', variable=self.datos,
                        value='total').grid(column=3, row=1, sticky=(tk.W), padx=5, pady=1)
        tabla = ttk.Labelframe(self, text='Models table', name='t_comp')
        tabla.grid(column=0, row=2, sticky=(tk.E, tk.W), padx=5, pady=5)
        self.tree = ttk.Treeview(tabla, height=12)
        self.tree['show'] = 'headings'
        sb = ttk.Scrollbar(tabla, orient="vertical", command=self.tree.yview)
        sb.grid(row=1, column=1, sticky="EW", pady=5)
        self.tree.configure(yscrollcommand=sb.set)
        self.tree["columns"] = ("1", "2", "3", "4")
        self.tree.column("1", width=70)
        self.tree.column("2", width=50)
        self.tree.column("3", width=350)
        self.tree.column("4", width=350)
        self.tree.heading("1", text="Nombre")
        self.tree.heading("2", text="size")
        self.tree.heading("3", text="variables")
        self.tree.heading("4", text="coeficientes")
        item = self.tree.insert("", "end", values=("",))
        self.tree.grid(row=1, column=0, padx=5, pady=5)
        frame = ttk.Labelframe(self, text='Select contribution types to visualise', name='contributions')
        frame.grid(column=0, row=3, columnspan=7, sticky=(tk.E, tk.W), padx=5, pady=5)
        self.csv_path.trace('w', self.__csv_path_changed)
        # select properties
        overall_state = True

        ttk.Radiobutton(frame, text='Overall', variable=self.contr_name,
                        value='total').grid(column=0, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Standard Distance', variable=self.contr_name,
                        value='Std').grid(column=1, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Dipole moments', variable=self.contr_name,
                        value='Dip').grid(column=2, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Dipole moments2', variable=self.contr_name,
                        value='Dip2').grid(column=3, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Hydrophobicity', variable=self.contr_name,
                        value='Hyd').grid(column=4, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Polar surface', variable=self.contr_name,
                        value='Pols').grid(column=5, row=0, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Molar refractivity', variable=self.contr_name,
                        value='Mol').grid(column=0, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Polarizability', variable=self.contr_name,
                        value='Pol').grid(column=1, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='van der waals radius', variable=self.contr_name,
                        value='Van').grid(column=2, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Gasteiger chargues', variable=self.contr_name,
                        value='Gas').grid(column=3, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Atomic weigth', variable=self.contr_name,
                        value='Ato').grid(column=4, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-R2', variable=self.contr_name,
                        value='Ab-R2').grid(column=5, row=1, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-pi2H', variable=self.contr_name,
                        value='Ab-pi2H').grid(column=0, row=2, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-sumA2H', variable=self.contr_name,
                        value='Ab-sumA2H').grid(column=1, row=2, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-sumB2H', variable=self.contr_name,
                        value='Ab-sumB2H').grid(column=2, row=2, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-sumB20', variable=self.contr_name,
                        value='Ab-sumB20').grid(column=3, row=2, sticky=(tk.W), padx=5, pady=1)
        ttk.Radiobutton(frame, text='Abraham-logL16', variable=self.contr_name,
                        value='Ab-logL16').grid(column=4, row=2, sticky=(tk.W), padx=5, pady=1)

        frame = ttk.Labelframe(self, text='Contribution view', name='vista')
        frame.grid(column=9, row=2, rowspan=2, sticky=(tk.EW, tk.N), padx=5, pady=5)
        self.original = Image.open(os.path.join(get_script_path(), 'No_mol.png'))
        resized = self.original.resize((450, 450), Image.ANTIALIAS)
        self.image = ImageTk.PhotoImage(resized)  # Keep a reference, prevent GC
        self.l2 = tk.Label(frame, image=self.image)
        self.l2.grid(row=0, column=0, columnspan=10)
        self.menos_btn = ttk.Button(frame, text='<', command=self.menos).grid(column=0, row=1, columnspan=2, padx=5, pady=5)
        self.label = tk.Label(frame, text="")
        self.label.grid(row=1, column=2, columnspan=6, padx=5, pady=5)
        self.mas_btn = ttk.Button(frame, text='>', command=self.mas).grid(column=8, row=1, columnspan=2, padx=5, pady=5)

        def cambiocheck(a):
            selectedItem = self.tree.selection()[0]
            self.n_mod = int(self.tree.item(selectedItem)['values'][0].split("_")[1])
            if self.tipo == 'bond':
                ruta = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                    "/Contr_" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "_" +
                                    self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "\\modelo_" + str(
                                        self.n_mod) + "\\" + str(self.contr_name.get()) + "_" + str(
                                        self.n_mol) + '.png')
            else:
                ruta = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                    "/Contr_ato" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "_" +
                                    self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "\\modelo_" + str(
                                        self.n_mod) + "\\" + str(self.contr_name.get()) + "_" + str(
                                        self.n_mol) + '.png')
            if os.path.exists(ruta):
                self.original = Image.open(ruta)
            else:
                self.original = Image.open(os.path.join(get_script_path(), 'No_mol.png'))
            resized = self.original.resize((400, 400), Image.ANTIALIAS)
            self.image = ImageTk.PhotoImage(resized)  # Keep a reference, prevent GC
            self.l2.config(image=self.image)
            self.label.config(text=self.names[self.n_mol-1])

        self.tree.bind("<<TreeviewSelect>>", cambiocheck)
        self.contr_name.trace('w', self.__contr_changed)
        #self.mas_btn.state(['disabled'])
        #self.menos_btn.state(['disabled'])

        # events
        #self.__save_fig_file_path.bind('<Return>', self.__run_plot)

        parent.add(self, text=tab_name)

    def mas(self):
        if self.tipo == 'bond':
            camino = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                  "/Contr_" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "_" + self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "\\modelo_" + str(self.n_mod) + "\\" + str(
                                      self.contr_name.get()) + "_" + str(self.n_mol+1) + '.png')
        else:
            camino = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                  "/Contr_ato" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "_" + self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[
                                      0] + "\\modelo_" + str(self.n_mod) + "\\" + str(
                                      self.contr_name.get()) + "_" + str(self.n_mol+1) + '.png')
        if os.path.exists(camino):
            self.n_mol += 1
            self.original = Image.open(camino)
        else:
            self.original = Image.open(os.path.join(get_script_path(), 'No_mol.png'))
        resized = self.original.resize((400, 400), Image.ANTIALIAS)
        self.image = ImageTk.PhotoImage(resized)  # Keep a reference, prevent GC
        self.l2.config(image=self.image)
        self.label.config(text=self.names[self.n_mol - 1])

    def menos(self):
        if (self.n_mol-1==0):
            self.n_mol = 1
        else:
            self.n_mol -= 1
        if self.tipo == 'bond':
            camino = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                  "/Contr_" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "_" +
                                  self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "\\modelo_" + str(
                                      self.n_mod) + "\\" + str(self.contr_name.get()) + "_" + str(
                                      self.n_mol) + '.png')
        else:
            camino = os.path.join(os.path.dirname(self.path.get()),"TOPSMODE" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0]+
                                  "/Contr_ato" + self.path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "_" +
                                  self.csv_path.get().rsplit('/', 1)[1].rsplit(".", 1)[0] + "\\modelo_" + str(
                                      self.n_mod) + "\\" + str(self.contr_name.get()) + "_" + str(
                                      self.n_mol) + '.png')
        if os.path.exists(camino):
            self.original = Image.open(camino)
        else:
            self.original = Image.open(os.path.join(get_script_path(), 'No_mol.png'))
        resized = self.original.resize((400, 400), Image.ANTIALIAS)
        self.image = ImageTk.PhotoImage(resized)  # Keep a reference, prevent GC
        self.l2.config(image=self.image)
        self.label.config(text=self.names[self.n_mol - 1])


def main():
    #platform: alt, clam, classic, default
    #Windows: vista, winnative, xpnative

    root = tk.Tk()
    #style = ttk.Style(root)
    root.tk.call('source', 'C:/DESCRIPTORES/Modeslab_gui/Azure-ttk-theme-main/azure.tcl')
    root.tk.call("set_theme", "dark")
    #root.tk.call('source', 'C:/script_python/Azure-ttk-theme-main/azure.tcl')
    #style.theme_use('azure')
    #style.theme_use('aqua')
    root.title("TOPSMODE - TOPological Sub-Structural MOlecular DEsign interpretation of QSAR models")
    root.geometry("1600x900")
    #root.iconbitmap('hola.ico')
    content = ttk.Frame(root)
    sdf_frame = Sdfreader(content, 'Molecules input', 0, 0, True)
    #content.config(bg="lightblue")
    lbl_copyright = ttk.Label(content, text='(c) Alfonso Pérez-Garrido 2023')

    tabs = ttk.Notebook(content)
    tab_1 = Tab_1(tabs, 'Calculate Descriptors', path=sdf_frame.sdf_path, nombres=sdf_frame.compound_names, act_field=sdf_frame.activity_field_name)
    tab_2 = Tab_2(tabs, 'Plot contributions', path=sdf_frame.sdf_path, nombres=sdf_frame.compound_names)
    content.grid(column=0, row=0, sticky=(tk.W,  tk.E))
    sdf_frame.grid(column=0, row=1, sticky=(tk.W, tk.N, tk.E))
    tabs.grid(column=0, row=2, sticky=(tk.W, tk.E))
    lbl_copyright.grid(column=0, row=3, sticky=(tk.W, tk.E, tk.S))

    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    content.columnconfigure(0, weight=1)
    content.rowconfigure(0, weight=1)

    root.mainloop()

if __name__ == '__main__':
    main()