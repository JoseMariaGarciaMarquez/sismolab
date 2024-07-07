import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sismolab_app import SeismogramApp

def test_load_data():
    app = SeismogramApp()
    app.load_data()
    assert app.st is not None
    assert len(app.st) > 0