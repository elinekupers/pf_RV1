import time

def time_string():
    """
    Returns a time string in local time
    Year_Month_Day_Hour:Minute:Second
    """
    return time.strftime("%Y_%m_%d__%H_%M_%S", time.localtime())