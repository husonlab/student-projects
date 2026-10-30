"""MAGrix: query large MAG x gene family matrices without loading them into memory.

    from magrix import PanGenomeMatrices
    m = PanGenomeMatrices("output/my_dataset")
    families, copies = m.paralog.row(m.dicts.mag_id("Mag568575"))
"""
from .api import Container, Dictionaries, MembersView, PanGenomeMatrices, PresenceView, SparseMatrix

__version__ = "0.1.0"

__all__ = [
    "PanGenomeMatrices",
    "SparseMatrix",
    "PresenceView",
    "MembersView",
    "Dictionaries",
    "Container",
]
