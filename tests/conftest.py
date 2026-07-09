from pathlib import Path

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture
def nopal_df():
    return pd.read_csv(REPO_ROOT / "nopal.csv")


@pytest.fixture
def nopal_nutrient_cols():
    return ["N", "P", "K", "Ca", "Mg"]


@pytest.fixture
def maiz_df():
    return pd.read_csv(REPO_ROOT / "maiz.csv")


@pytest.fixture
def maiz_nutrient_cols(maiz_df):
    return [c for c in maiz_df.select_dtypes("number").columns if c != "yield"]
