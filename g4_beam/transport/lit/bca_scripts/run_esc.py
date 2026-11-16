import subprocess, os, sys
base = "."  # directory that holds your RustBCA checkout
exe = base + "/RustBCA/target/release/RustBCA"
run = base + "/bca_runs"
os.chdir(run)
targets = {
 "Ni":  dict(n=[9.14e10], Es=[4.44], Ec=[3.0], Eb=[0.0], Z=[28], m=[58.69]),
 "W":  dict(n=[6.31e10], Es=[8.79], Ec=[3.0], Eb=[0.0], Z=[74], m=[183.84]),
 "Pd":  dict(n=[6.80e10], Es=[3.89], Ec=[3.0], Eb=[0.0], Z=[46], m=[106.42]),
 "Ti":  dict(n=[5.67e10], Es=[4.84], Ec=[3.5], Eb=[0.0], Z=[22], m=[47.867]),
 #"TiO2":dict(n=[3.19e10, 6.38e10], Es=[4.84, 2.58], Ec=[3.5, 2.0], Eb=[3.0, 2.58], Z=[22, 8], m=[47.867, 15.9994]),
}
energies = [500, 1000, 6500, 13000]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
tmpl = """[options]
name = "{name}"
track_trajectories = false
track_recoils = true
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
electronic_stopping_correction_factor = 1.53

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
print("target E_eV R_N Y_total Y_by_mass")
for t, p in targets.items():
    for E in energies:
        name = f"ESC_D_{t}_{E}_"
        txt = tmpl.format(name=name, N=N, E=float(E), n=p["n"], Eb=p["Eb"], Es=p["Es"], Ec=p["Ec"], Z=p["Z"], m=p["m"], ii=[0]*len(p["Z"]))
        fn = name + ".toml"
        open(fn, "w").write(txt)
        subprocess.run([exe, "0D", fn], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        rows = [l.strip().rstrip(",").split(",") for l in open(name + "summary.output").read().splitlines()[1:] if l.strip()]
        RN = 0; Ys = {}
        for r in rows:
            mass, refl, sput, dep = float(r[0]), int(r[1]), int(r[2]), int(r[3])
            if abs(mass - 2.014) < 0.01: RN = refl / N
            else: Ys[round(mass,1)] = sput / N
        print(t, E, f"{RN:.3f}", f"{sum(Ys.values()):.4f}", {k: round(v,4) for k,v in Ys.items()}, flush=True)
