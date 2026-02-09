import luteastress_native
import time
import numpy as np

print("Running NATIVE high-fidelity mechanistic trial...")
# We need to see the raw ovulation times to debug OSI=1.0
results = luteastress_native.run_mechanistic_experiment(
    n_individuals=3,
    n_days=100,
    n_trials=1,
    seed=42
)

# Wait, the native code doesn't return the raw ovulation times, just OSI.
# I need to modify the native code to return them for debugging, 
# or just add some print statements in Rust.

print(f"Treatment OSI: {results['treatment_osi']}")
print(f"Control OSI: {results['control_osi']}")
