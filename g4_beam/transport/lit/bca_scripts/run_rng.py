import subprocess, os, sys
base = "."  # directory that holds your RustBCA checkout
exe = base + "/RustBCA/target/release/RustBCA"
run = base + "/bca_runs"
os.chdir(run)
targets = {
 "Pd":  dict(n=[6.80e10], Es=[3.89], Ec=[3.0], Eb=[0.0], Z=[46], m=[106.42]),
 "Ti":  dict(n=[5.67e10], Es=[4.84], Ec=[3.5], Eb=[0.0], Z=[22], m=[47.867]),
 #"TiO2":dict(n=[3.19e10, 6.38e10], Es=[4.84, 2.58], Ec=[3.5, 2.0], Eb=[3.0, 2.58], Z=[22, 8], m=[47.867, 15.9994]),
}
energies = [250, 500, 1000, 2500, 5000, 6500, 13000]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
tmpl = """[options]
name = "{name}"
track_trajectories = false
track_recoils = false
track_recoil_trajectories = false
write_buffer_size = 8000
weak_collision_order = 3
suppress_deep_recoils = false
high_energy_free_flight_paths = false
electronic_stopping_mode = "LOW_ENERGY_NONLOCAL"
mean_free_path_model = "LIQUID"
interaction_potential = [["KR_C"]]
scattering_integral = [["MENDENHALL_WELLER"]]
root_finder = [[{{"NEWTON"={{max_iterations=100, tolerance=1E-6}}}}]]
num_threads = 8
num_chunks = 50
track_displacements = false
track_energy_losses = false
energy_min = 0.0
energy_max = 20000.0
energy_num = 10
angle_min = 0.0
angle_max = 90.0
angle_num = 10
x_min = 0.0
y_min = -0.1
z_min = -0.1
x_max = 0.1
y_max = 0.1
z_max = 0.1
x_num = 10
y_num = 10
z_num = 10

[particle_parameters]
length_unit = "1e-6"
energy_unit = "EV"
mass_unit = "AMU"
N = [ {N},]
m = [ 2.014,]
Z = [ 1,]
E = [ {E},]
Ec = [ 0.1,]
Es = [ 0.0,]
interaction_index = [0]
pos = [ [ 0.0, 0.0, 0.0,],]
dir = [ [ 0.9999, 0.0001, 0.0,],]
particle_input_filename = ""

[geometry_input]
length_unit = "MICRON"
densities = {n}
electronic_stopping_correction_factor = {esc}

[material_parameters]
energy_unit = "EV"
mass_unit = "AMU"
Eb = {Eb}
Es = {Es}
Ec = {Ec}
Z = {Z}
m = {m}
interaction_index = {ii}
surface_binding_model = {{"PLANAR"={{calculation="TARGET"}}}}
bulk_binding_model = "AVERAGE"
"""
print("target E_eV esc R_N mean_depth std median")
for t, p in targets.items():
    for E in energies:
     for esc in (1.0, 1.53):
        name = f"RNG_{t}_{E}_{esc}_"
        txt = tmpl.format(esc=esc, name=name, N=N, E=float(E), n=p["n"], Eb=p["Eb"], Es=p["Es"], Ec=p["Ec"], Z=p["Z"], m=p["m"], ii=[0]*len(p["Z"]))
        fn = name + ".toml"
        open(fn, "w").write(txt)
        subprocess.run([exe, "0D", fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        import numpy as np
        d = np.loadtxt(name + "deposited.output", delimiter=",", ndmin=2)
        x = d[:, 2] * 1000.0  # micron -> nm
        RN = 1 - len(x) / N
        print(t, E, esc, f"R_N={RN:.3f}", f"Rp={x.mean():.2f} nm", f"dRp={x.std():.2f} nm", f"median={np.median(x):.2f} nm", flush=True)
