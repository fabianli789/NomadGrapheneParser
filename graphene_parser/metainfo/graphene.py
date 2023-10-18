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

from nomad.metainfo import MSection, Section, SubSection, Quantity, Package
from nomad.datamodel.metainfo import simulation
from nomad.datamodel import results, optimade
#from nomad.datamodel.metainfo.workflow import Workflow


m_package = Package()


class Dimensions_Graphene(MSection):
    m_def = Section(validate=False)
    x = Quantity(type=float, description='x-dimension of 2d lattice in nm')
    y = Quantity(type=float, description='y-dimension of 2d lattice in nm')
    mean_radius = Quantity(type=float, description='Mean Radius of Graphene sample in nm')
    ra = Quantity(type=float, description='Roughness RA in nm')
    rq = Quantity(type=float, description='Roughness RQ in nm')
    hydro_edge = Quantity(type=float, description='Share of Hydrogenated edges in %')
    defects = Quantity(type=float, description='Share of defects in %')
    
class Concentrations_Graphene(MSection):
    m_def = Section(validate=False)
    name = Quantity(type=str, description='name of molecule. "e" means already attached to graphene edge.')
    concentration = Quantity(type=float, shape = ['*'], description='concentration of molecule at specific time, same length as "time"-array at run.calculation.concentration_time')


class ChemReactions_Graphene(MSection):
    name = Quantity(type=str, description = 'name of chem. reaction')
    barrier = Quantity(type=float, shape=[], description='energetic barrier in eV')
    
    occurences = Quantity(type=int, shape=[], description ='number of occurences of this chem. reaction')
    residence_time = Quantity(type=float, shape=[], description =  'time of each chem. reaction')
    

class Calculation(simulation.calculation.Calculation):
    m_def = Section(validate=False, extends_base_section=True)    
    dimensional_properties = SubSection(sub_section=Dimensions_Graphene.m_def, repeats=False)
    
    chem_reactions = SubSection(sub_section=ChemReactions_Graphene.m_def, repeats=True)
    concentrations = SubSection(sub_section=Concentrations_Graphene.m_def, repeats=True)
    volume_fraction = Quantity(type=float, description='volume if SEI got pressed together')
    porosity = Quantity(type=float, description='share of porous volume wrt to total SEI volume')
    concentration_time = Quantity(type=float, shape=['*'], description='time evolution for the concentration of molecules, same length as "concentration"-array under run.calculation.concentrations')
    pressure_CH4 = Quantity(type=float, description='Pressure of CH4 in Torr')
    pressure_H2 = Quantity(type=float, description='Pressure of H2 in Torr')
    bonds  = Quantity(type=int, shape=['*'], description = 'index of paired atom. -1 means unavailable. Same order as "cartesian_site_coordinates" - array.')
    flag = Quantity(type=int, shape=['*'], description = '# of bonds in graphene flake. Same order as "cartesian_site_coordinates" - array.')
    mean_radius_growth = Quantity(type=float, shape=['*'], description = 'change of mean radius  in nm over time. See mean_radius_growth_time for time steps.')
    mean_radius_growth_time = Quantity(type=float, shape=['*'], description = 'time steps for mean radius growth. Same order of array as "mean_radius_growth".')
class Run(simulation.run.Run):
    m_def = Section(validate=False, extends_base_section=True)

m_package.__init_metainfo__()
