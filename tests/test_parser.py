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

import pytest
import logging

from nomad.datamodel import EntryArchive

from graphene_parser import GrapheneParser


@pytest.fixture
def parser():
    return GrapheneParser()


def test_example(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Calculation.in', archive, logging)

    run = archive.run[0]
    calc = run.calculation[0]
    mol = calc.mols[0]
    assert mol.name == "A.pqr"
