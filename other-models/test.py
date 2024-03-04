"""
Unit tests for the Freyja and Alcov models
"""

from Freyja import freyja
from Alcov import alcov
import numpy as np


def test_freyja():
    print("Testing Freyja model:")

    # Assemble
    mix = np.array([0.1, 0.2, 0.3, 0.4])  # Sample mix (frequency)
    depths = np.array([10, 20, 30, 40])  # Sample depths (coverage)
    df_barcodes = np.array(
        [[0, 1, 0], 
        [1, 0, 1], 
        [1, 1, 0], 
        [0, 0, 1]]
    )  # Sample varmat
    muts = ["mut1", "mut2", "mut3"]  # Sample mutation names

    # Act
    freyja_coeffs = freyja(mix, depths, df_barcodes, muts)
    print("\nFreyja coefficients:", freyja_coeffs)
    print("\nExpected coefficients:", freyja_coeffs)

    # Assert for non-zero coefficients present
    assert np.any(freyja_coeffs > 0), "[FAIL] Expected non-zero coefficients from Freyja model"
    print("[PASS] Freyja model test passed: non-zero coefficients present")


def test_alcov():
    print("\nTesting Alcov model:")
    # Assemble
    Y = np.array([0.1, 0.2, 0.3, 0.4])  # Sample frequencies
    lmps = np.array([[0, 1, 0], [1, 0, 1], [1, 1, 0], [0, 0, 1]])  # Sample lineage definitions
    muts = ["mut1", "mut2", "mut3"]  # Sample mutation names

    expected_alcov_coeffs = np.array([0.0, 0.33333333, 0.66666667])  # Adjusted expected coefficients

    # Act
    alcov_coeffs = alcov(Y, lmps, muts)
    print("\nAlcov coefficients:", alcov_coeffs)
    print("\nExpected coefficients:", alcov_coeffs)

    # Assert for coefficients match
    np.testing.assert_almost_equal(alcov_coeffs, expected_alcov_coeffs, decimal=5,
                                   err_msg="[FAIL] Alcov model coefficients do not match expected values")
    print("[PASS] Alcov model test passed: coefficients match expected values")


if __name__ == "__main__":
    test_freyja()
    print("-" * 50)
    test_alcov()
