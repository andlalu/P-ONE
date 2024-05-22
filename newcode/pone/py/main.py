# Sandbox to test the particle filtering

from ccf import cf_coeffs_heston
import ponecpp

def test():
    testRandInstance = ponecpp.Rand()
    print(testRandInstance.rand_int())

if __name__ == '__main__':
    test()

