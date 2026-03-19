from pathlib import Path


CASE_ORDER = (
    "traction",
    "water_loading",
    "des3d_traction",
    "des3d_water_loading",
)


CASE_CONFIG = {
    "traction": {
        "cfg": Path("1d-consolidation-des2d_traction.cfg"),
        "model": "terzaghi_traction",
        "label": "Traction loading",
        "dim": 2,
        "skip_first_output": 1,
        "run_subdir": "traction",
    },
    "water_loading": {
        "cfg": Path("1d-consolidation-des2d_water_loading.cfg"),
        "model": "terzaghi_water_loading",
        "label": "Water loading",
        "dim": 2,
        "skip_first_output": 0,
        "run_subdir": "water_loading",
    },
    "des3d_traction": {
        "cfg": Path("1d-consolidation-des3d_traction.cfg"),
        "model": "terzaghi",
        "label": "Traction loading (3D)",
        "dim": 3,
        "skip_first_output": 1,
        "run_subdir": "des3d_traction",
    },
    "des3d_water_loading": {
        "cfg": Path("1d-consolidation-des3d_water_loading.cfg"),
        "model": "terzaghi",
        "label": "Water loading (3D)",
        "dim": 3,
        "skip_first_output": 0,
        "run_subdir": "des3d_water_loading",
    },
}


def case_names():
    return list(CASE_ORDER)


def case_config(case_name):
    return CASE_CONFIG[case_name]
