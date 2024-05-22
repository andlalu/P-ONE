from ctypes import *
import os

shared_library_dir = "/Users/andrei/P-ONE/newcode/pone/cpp/build/debug"
cdll.LoadLibrary(os.path.join(shared_library_dir, "pone.cpython-311-darwin.so"))



