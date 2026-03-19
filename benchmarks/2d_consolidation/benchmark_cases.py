from pathlib import Path


CASE_ORDER = (
    "mandel",
)


CASE_CONFIG = {
    "mandel": {
        "cfg": Path("2d-consolidation.cfg"),
        "model": "mandel",
        "label": "Mandel consolidation",
        "dim": 2,
        "run_subdir": "mandel",
    },
}


def case_names():
    return list(CASE_ORDER)


def case_config(case_name):
    return CASE_CONFIG[case_name]
