#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import re
import datetime
import numpy as np
from pathlib import Path

from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser
from nomad.units import ureg as units
from nomad.datamodel.metainfo.simulation.run import Run, Program
from nomad.datamodel.metainfo.simulation.system import System
from nomad.datamodel.metainfo.simulation.calculation import Calculation, Energy, EnergyEntry
#from nomad.datamodel.metainfo.workflow import Workflow
from nomad.datamodel.results import Results, Properties, Structure
from nomad.parsing.file_parser import UnstructuredTextFileParser, Quantity
from nomad.datamodel.optimade import Species
from . import metainfo  # pylint: disable=unused-import
from .metainfo.graphene import Dimensions_Graphene, ChemReactions_Graphene, Concentrations_Graphene


def DetailedParser(filepath, archive):
    
    if os.path.exists(str(filepath.parent) + r'/status_graphene.csv'):
        with open(str(filepath.parent) + r'/status_graphene.csv') as status_file:
            time_run = archive.m_setdefault("run.time_run")
            time_run.cpu1_start = 0
            calc = archive.m_setdefault("run.calculation")
        
            for i, line in enumerate(status_file):
                line = line.strip("\n")
                parts = line.split(",")
                if parts[0] == None:
                    continue
                if re.search(r"cpu", parts[0].lower()):
                    time_run.cpu1_end = parts[1]
                if re.search(r'kmc time', parts[0].lower()):
                    calc.time = parts[1]
                if re.search(r'step', parts[0].lower()):
                    calc.step = int(float(parts[1]))
    

    if os.path.exists(str(filepath.parent) + r'/occurrence_res_graphene.csv'):
        with open(str(filepath.parent) + r"/occurrence_res_graphene.csv") as occurrence_file:
            occurence_array = []
            residence_time_array = []
            for i, line in enumerate(occurrence_file):
                if re.search("Oc", line):
                    continue
                parts = line.split(",")
                occurence_array.append(int(parts[1]))
                if len(parts)>=3:
                    residence_time_array.append(float(parts[2]))
                else:
                    pass
    if os.path.exists(str(filepath.parent) + r'/input_graphene.yml'):    
        with open(str(filepath.parent) + r'/input_graphene.yml') as file:
            dim = calc.m_create(Dimensions_Graphene)
            j = 0
            for i, line in enumerate(file):
                parts  = line.split(": ")
               
                if re.search(r'x_dim', parts[0].lower()):
                    dim.x = float(parts[1])
                elif re.search(r'y_dim', parts[0].lower()):
                    dim.y = float(parts[1])
                elif re.search(r'p_ch4', parts[0].lower()):
                    calc.pressure_CH4 = float(parts[1])
                elif re.search(r'p_h2', parts[0].lower()):
                    calc.pressure_H2 = float(parts[1])        
                elif re.search(r'\-', parts[0]):
                    chem_reactions = calc.m_create(ChemReactions_Graphene)
                    parts[0] = parts[0].lstrip('- ').rstrip(' ')
                    chem_reactions.name = parts[0]
                    chem_reactions.barrier = float(parts[1])
                    
                    chem_reactions.occurences = occurence_array[i-9]
                    if len(residence_time_array) > 0:
                        chem_reactions.residence_time = residence_time_array[i-9]
                    else:
                        pass
    if os.path.exists(str(filepath.parent) + r'/last_step_graphene.csv'):
        with open(str(filepath.parent) + r'/last_step_graphene.csv') as last_step_file:
            species_array = []
            for j, x in enumerate(last_step_file):
                pass
            bonds = []
            flag = []
            coord_x_final = []
            coord_y_final = []
            coordinates_final =  np.zeros((j,3))
            
            last_step_file.seek(0)    
            for i, line in enumerate(last_step_file):
                parts = line.strip("\n").split(",")
                if re.search(r'species', line):
                    continue
                if re.search(r'\d', parts[0]):
                    coord_x_final.append(float(parts[0].strip('"').replace('[', '')))
                else:
                    coord_x_final.append(-1)
                if re.search(r'\d', parts[1]):
                    coord_y_final.append(float(parts[1].strip('"').replace(']', '')))
                else:
                    coord_y_final.append(-1)

                species_array.append(parts[2])
 
                if re.search(r'\d', parts[3]):
                    bonds.append(int(parts[3].strip('"').replace('[', '').replace(']', '')))  
                else:
                    bonds.append(-1)
                if re.search(r'\d', parts[4]):
                   flag.append(int(parts[4].strip('"').replace('[', '').replace(']', '')))  
                else:
                   flag.append(-1)
            
            
            calc.bonds = np.array(bonds)
            calc.flag = np.array(flag)
            coordinates_final[:, 0] = np.array(coord_x_final)
            coordinates_final[:, 1] = np.array(coord_y_final)
            
            calc.species_coordinates = coordinates_final
            calc.species = species_array

                  
    if os.path.exists(str(filepath.parent) + r'/concentration_graphene.csv'):
        with open(str(filepath.parent) + r'/concentration_graphene.csv') as conc_file:
            first_line_parts = conc_file.readline().strip("\n").split(",")
            for x, bla in enumerate(conc_file):
                rows = x
            conc_array = np.zeros((rows+1, len(first_line_parts)-1))
            time_array = []
            
    
            conc_file.seek(0) 
            for j, line in enumerate(conc_file):
            
                if re.search(r'time', line.lower()):
                    continue
                
                parts = line.strip("\n").split(",")
                parts = [float(x) for x in parts]
                time_array.append(parts[0])
                parts = parts[1:]            
                conc_array[j-1] = parts
            
            calc.concentration_time = time_array

            for i in range(len(first_line_parts)-1):
                conc = calc.m_create(Concentrations_Graphene)   
                conc.name = first_line_parts[i+1]
                conc.concentration = conc_array[:,i]
    if os.path.exists(str(filepath.parent) + r'/growth_graphene.csv'):
        with open(str(filepath.parent) + r'/growth_graphene.csv') as growth_file:
            mean_radius_growth = []
            mean_radius_growth_time = []
            for i, line in enumerate(growth_file):
                if re.search(r'time', line.lower()):
                    continue
                parts = line.strip('"').split(',')
                mean_radius_growth.append(float(parts[1]))
                mean_radius_growth_time.append(float(parts[2]))        
            calc.mean_radius_growth = np.array(mean_radius_growth)
            calc.mean_radius_growth_time = np.array(mean_radius_growth_time)

    if os.path.exists(str(filepath.parent) + r'/properties_graphene.csv'):
        with open(str(filepath.parent) + r'/properties_graphene.csv') as prop_file:
            dim = calc.m_create(Dimensions_Graphene)
            for i, line in enumerate(prop_file):
                line = line.strip("\n")
                parts = line.split(",")
                if parts[0] == None:
                    continue
                if re.search(r'mean', parts[0].lower()):
                    dim.mean_radius = float(parts[1])    
                if re.search(r'ra', parts[0].lower()):
                    dim.ra = float(parts[1])
                if re.search(r'rq', parts[0].lower()):
                    dim.rq = float(parts[1])
                if re.search(r'hy', parts[0].lower()):
                    dim.hydro_edge = float(parts[1])
                if re.search(r'defe', parts[0].lower()):
                    dim.defects = float(parts[1])    
class GrapheneParser():

    def parse(self, filepath, archive, logger):
        sec_run = archive.m_create(Run)
        sec_program = archive.m_setdefault('run.program')
        sec_program.name = 'Meysam Testing cron job'
        sec_method = archive.m_setdefault('results.method')
        sec_method.method_name = 'unavailable'
        mainfile = Path(filepath)
        DetailedParser(mainfile, archive)
