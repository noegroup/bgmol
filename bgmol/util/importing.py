

__all__ = ["import_openmm"]


def import_openmm():
    try:
        import openmm as mm
        import openmm.unit as unit
        import openmm.app as app
    except ImportError:
        # fall back to older version < 7.6
        from simtk import openmm as mm
        from simtk import unit as unit
        from simtk.openmm import app as app
    return mm, unit, app