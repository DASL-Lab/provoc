import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import MinMaxScaler


def alcov(Y, lmps, muts):
    """
    A simplified Alcov model for linear regression without intercepts,
    ensuring positivity of the coefficients

    Parameters:
    - Y: Frequencies (numpy array or pandas series)
    - lmps: Lineage definitions, analogous to varmat (numpy array or pandas DataFrame)
    - muts: Mutation names (list or numpy array)

    Returns:
    - Coefficients of the linear regression model
    """

    # Ensure Y, lmps, and muts are properly aligned, crucial as Alcov ensures mutations match up
    if isinstance(lmps, pd.DataFrame):
        lmps = lmps[muts].to_numpy()
    else:
        # Assuming lmps is already filtered to match 'muts'
        pass

    # Scale Y and lmps for better regression performance
    scaler = MinMaxScaler()
    Y_scaled = scaler.fit_transform(Y.reshape(-1, 1)).flatten()
    lmps_scaled = scaler.fit_transform(lmps)

    # Linear Regression without intercept
    model = LinearRegression(fit_intercept=False, positive=True)
    model.fit(lmps_scaled, Y_scaled)

    # Extract and return the coefficients
    return model.coef_
