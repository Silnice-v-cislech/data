from datetime import datetime

from workflow.scripts.utils import zero_monthly_time_series


def test_zero_monthly_time_series():
    assert zero_monthly_time_series(datetime(2012, 12, 1), datetime(2014, 3, 30)) == {
        (2012, 12): 0,
        (2013, 1): 0,
        (2013, 2): 0,
        (2013, 3): 0,
        (2013, 4): 0,
        (2013, 5): 0,
        (2013, 6): 0,
        (2013, 7): 0,
        (2013, 8): 0,
        (2013, 9): 0,
        (2013, 10): 0,
        (2013, 11): 0,
        (2013, 12): 0,
        (2014, 1): 0,
        (2014, 2): 0,
        (2014, 3): 0,
    }


def test_zero_monthly_time_series_same_year():
    assert zero_monthly_time_series(datetime(2013, 3, 20), datetime(2013, 5, 15)) == {
        (2013, 3): 0,
        (2013, 4): 0,
        (2013, 5): 0,
    }
