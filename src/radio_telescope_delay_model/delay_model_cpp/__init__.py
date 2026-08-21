"""CALC11 delay model extension (C++ driver + reference Fortran core)."""

from __future__ import annotations


def calc_available() -> bool:
    """Whether the C++/Fortran extension was built and can be imported."""
    try:
        from radio_telescope_delay_model.delay_model_cpp import (  # noqa: F401
            _delay_model_ext,
        )
    except ImportError:
        return False
    return True


def almacalc(*args, **kwargs):
    """In-place CALC11 delay computation; see ``_delay_model_ext.almacalc``."""
    from radio_telescope_delay_model.delay_model_cpp import _delay_model_ext

    return _delay_model_ext.almacalc(*args, **kwargs)
