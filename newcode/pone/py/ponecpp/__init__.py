import os
import sys
import importlib.util

# Path to the compiled shared library
current_dir = os.path.dirname(__file__)
parent_dir = os.path.dirname(current_dir)
grandparent_dir = os.path.dirname(parent_dir)
lib_name = 'pone.cpython-311-darwin.so'  # Make sure this matches your actual file name
lib_path = os.path.join(grandparent_dir, 'cpp/build/debug/', lib_name)

# Load the shared library
spec = importlib.util.spec_from_file_location('pone', lib_path)
pone = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pone)

# Import the entire module
sys.modules['ponecpp.pone'] = pone



