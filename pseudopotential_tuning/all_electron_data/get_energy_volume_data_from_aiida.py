# Data: https://archive.materialscloud.org/record/2023.81

import numpy as np
import aiida
from aiida.orm import load_group
from aiida.common import LinkType
import json

aiida.load_profile()


target_element = 'Li'

target_configuration = 'X/BCC'
WORKFLOWS_GROUP_LABEL = f'acwf-verification/unaries-set2/workflows/wien2k'

# target_configuration = 'XO'
# WORKFLOWS_GROUP_LABEL = f'acwf-verification/oxides-verification-PBE-v1/workflows/wien2k'

group = load_group(WORKFLOWS_GROUP_LABEL)

for node in group.nodes:
    structure = node.inputs.structure

    # Check if structure is Li-BCC
    if structure.base.extras.all['element'] != target_element or structure.base.extras.all['configuration'] != target_configuration:
        continue

    # Collect all volumes and energies for this system
    volumes = []
    energies = []
    # Filter successful workflows
    if node.process_state.value == 'finished' and node.exit_status == 0:
        # Get all output links of the workflow, of type return
        outputs = node.base.links.get_outgoing(link_type=LinkType.RETURN).nested()
        # Loop over all output structures, get the volume and the corresponding
        # energy
        for index, sub_structure in sorted(outputs['structures'].items()):
            volumes.append(sub_structure.get_cell_volume())
            energies.append(outputs['total_energies'][index].value)

    # Sort (V, E) pairs
    energies = [e for _, e in sorted(zip(volumes, energies))]
    volumes = sorted(volumes)
    # # print volume and energy
    # for V, E in zip(volumes, energies):
    #     print(f"{V} {E}")

    print(json.dumps({'volumes_Ang3': volumes, 'energies_eV': energies}, indent=4))

    print()
    break
