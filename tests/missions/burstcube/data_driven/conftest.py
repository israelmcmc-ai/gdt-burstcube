# Copyright 2024-2025 by the BurstCube Team.
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
"""Shared helpers for the data-driven tests. These need the real files
listed in ``src/gdt/data/burstcube.urls``, downloaded with
``gdt-data download burstcube`` into ``gdt.core.data_path / 'burstcube'``;
any test whose file is missing skips cleanly rather than failing.
"""
from pathlib import Path

import pytest

from gdt.core import data_path

BURSTCUBE_DATA_DIR = Path(data_path) / 'burstcube'


def real_file(basename: str) -> Path:
    """Resolve a real downloaded data file by basename, skipping the test
    if it is not present.

    Args:
        basename (str): The file's basename, as listed in
            ``src/gdt/data/burstcube.urls``

    Returns:
        (Path)
    """
    path = BURSTCUBE_DATA_DIR / basename
    if not path.is_file():
        pytest.skip(f'{basename} not found in {BURSTCUBE_DATA_DIR} -- run '
                   f'`gdt-data download burstcube` to fetch the data-driven '
                   f'test files.')
    return path
