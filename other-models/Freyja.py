from sklearn.linear_model import Lasso
import numpy as np

def freyja(mix, depths, df_barcodes, muts, eps=1e-4):
    """
    Simplified Freyja model using Lasso regression.

    Parameters:
    - mix: Array of frequency (count divided by coverage) for each sample.
    - depths: Array of coverage for each sample.
    - df_barcodes: 2D array representing varmat, with rows as samples and columns as mutations.
    - muts: List of mutation names.
    - eps: Regularization strength for the Lasso regression.

    Returns:
    - Coefficients from the Lasso regression, representing the estimated proportions of variants.
    """

    # Adjust the importance of mutations based on coverage
    depth_adjustment = np.log(depths + 1) / np.max(np.log(depths + 1))
    adjusted_mix = mix * depth_adjustment

    # Apply depth adjustment to df_barcodes
    adjusted_barcodes = df_barcodes * depth_adjustment[:, np.newaxis]

    # Initialize and fit the Lasso model without intercepts and ensuring positivity
    lasso = Lasso(alpha=eps, fit_intercept=False, positive=True, max_iter=10000)
    lasso.fit(adjusted_barcodes, adjusted_mix)

    return lasso.coef_
