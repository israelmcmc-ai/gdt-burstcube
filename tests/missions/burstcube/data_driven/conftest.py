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
