"""virtual-macrophage — a single-cell macrophage agent.

v0.1: deterministic core, perceive/decide/act loop, literature-grounded rules.
"""
from .macrophage import Macrophage
from .environment.world import World

__version__ = "0.1.0"
__all__ = ["Macrophage", "World"]
