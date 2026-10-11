# model/parameters.py
# Migration costs κ^a, calibrated in config.py (model/calibration.py) so that the
# symmetric baseline reproduces the staying probabilities implied by the data.
from config import κl_cal, κh_cal

κl, κh = κl_cal, κh_cal
